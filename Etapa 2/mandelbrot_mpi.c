/* Comandos para compilar e executar (é necessário entrar na pasta da Etapa 2)
 *
 * Compilar (no Windows):
 * gcc mandelbrot_mpi.c -O3 -std=c11 -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include" -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi -lm -o mandelbrot_mpi.exe
 *
 * Executar:
 * mpiexec -n 4 .\mandelbrot_mpi.exe blocos
 * mpiexec -n 4 .\mandelbrot_mpi.exe ciclica
 */
#include <mpi.h>
#include <limits.h>
#include <errno.h>
#define BASE_FALHA() MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#define OUTPUT_DIR "saida_mpi"

typedef struct {
  unsigned char r;
  unsigned char g;
  unsigned char b;
} RGB;

// Cria a pasta de saída
void cria_pasta_saida() {
#ifdef _WIN32
  _mkdir(OUTPUT_DIR);
#else
  mkdir(OUTPUT_DIR, 0755);
#endif
}

// Calcula a matriz sequencial para validação
int32_t *mandelbrot(double re_min, double re_max, double im_min, double im_max, int width, int height, int max_iter) {

  size_t total = (size_t)width * (size_t)height;
  int32_t *count = malloc(total * sizeof(int32_t));

  if (count == NULL) {
    fprintf(stderr, "Erro: nao foi possivel alocar memoria.\n");
    BASE_FALHA();
  }

  for (int py = 0; py < height; py++) {
    double ci = im_min + ((double)py / (double)(height - 1)) * (im_max - im_min);

    for (int px = 0; px < width; px++) {
      double cr = re_min + ((double)px / (double)(width - 1)) * (re_max - re_min);
      double zr = 0.0;
      double zi = 0.0;
      int iteracoes = max_iter;

      for (int i = 0; i < max_iter; i++) {
        double novo_zr = zr * zr - zi * zi + cr;
        double novo_zi = 2.0 * zr * zi + ci;

        zr = novo_zr;
        zi = novo_zi;

        if (zr * zr + zi * zi > 4.0) {
          iteracoes = i;
          break;
        }
      }
      count[(size_t)py * width + px] = iteracoes;
    }
  }
  return count;
}

// Obtém a cor pela paleta selecionada
static RGB cor_paleta(double valor, int custo) {
  static const RGB malha[] = {
      {11, 29, 58}, {18, 56, 110}, {31, 111, 180}, {34, 182, 200}, {240, 180, 41}, {255, 243, 209}};
  static const RGB inferno[] = {
      {0, 0, 4}, {40, 11, 84}, {101, 21, 110}, {159, 42, 99}, {212, 72, 66}, {245, 125, 21}, {250, 193, 39}, {252, 255, 164}};
  const RGB *paleta = custo ? inferno : malha;
  int n = custo ? 8 : 6;
  if (valor < 0)
    valor = 0;
  if (valor > 1)
    valor = 1;
  double posicao = valor * (n - 1);
  int i = (int)posicao;
  if (i >= n - 1)
    return paleta[n - 1];

  double t = posicao - i;
  RGB a = paleta[i], b = paleta[i + 1];
  RGB cor = {
      (unsigned char)(a.r + (b.r - a.r) * t),
      (unsigned char)(a.g + (b.g - a.g) * t),
      (unsigned char)(a.b + (b.b - a.b) * t)};
  return cor;
}

// Grava a imagem em PPM
static void salva_ppm(const int32_t *imagem, int w, int h, int limite, int custo, const char *caminho) {
  FILE *f = fopen(caminho, "wb");
  if (!f)
    BASE_FALHA();
  fprintf(f, "P6\n%d %d\n255\n", w, h);
  double max_log = log1p((double)limite);

  for (int y = h - 1; y >= 0; y--) {
    for (int x = 0; x < w; x++) {
      int32_t iter = imagem[(size_t)y * w + x];
      RGB cor;
      if (custo)
        cor = cor_paleta((double)iter / limite, 1);
      else if (iter >= limite)
        cor = (RGB){(unsigned char)(0.043 * 255.0),
                    (unsigned char)(0.114 * 255.0),
                    (unsigned char)(0.227 * 255.0)};
      else
        cor = cor_paleta(log1p((double)iter) / max_log, 0);
      if (fwrite(&cor, sizeof cor, 1, f) != 1)
        BASE_FALHA();
    }
  }
  if (fclose(f))
    BASE_FALHA();
  printf("  Salvo: %s\n", caminho);
}

// Exibe o erro e encerra os processos
static void falha(const char *mensagem) {
  fprintf(stderr, "%s\n", mensagem);
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  exit(EXIT_FAILURE);
}

// Aloca memória e verifica falhas
static void *aloca(size_t n, size_t tamanho) {
  if (n > SIZE_MAX / tamanho)
    falha("Tamanho excede SIZE_MAX.");
  void *p = malloc((n ? n : 1) * tamanho);
  if (!p)
    falha("Memoria insuficiente.");
  return p;
}

// Converte e valida uma dimensão
static int inteiro(const char *s) {
  char *fim;
  errno = 0;
  long n = strtol(s, &fim, 10);
  if (errno || *fim || fim == s || n < 2 || n > INT_MAX)
    falha("Dimensoes devem ser inteiros entre 2 e INT_MAX.");
  return (int)n;
}

// Calcula as contagens de uma linha
static void linha(int32_t *destino, int py, int w, int h, int limite, double rmin, double rmax, double imin, double imax) {
  double ci = imin + ((double)py / (double)(h - 1)) * (imax - imin);
  for (int px = 0; px < w; px++) {
    double cr = rmin + ((double)px / (double)(w - 1)) * (rmax - rmin);
    double zr = 0.0, zi = 0.0;
    int iteracoes = limite;
    for (int i = 0; i < limite; i++) {
      double novo_zr = zr * zr - zi * zi + cr;
      double novo_zi = 2.0 * zr * zi + ci;
      zr = novo_zr;
      zi = novo_zi;
      if (zr * zr + zi * zi > 4.0) {
        iteracoes = i;
        break;
      }
    }
    destino[px] = iteracoes;
  }
}

// Exibe a configuração da execução
static void mostra_configuracao(int p, int w, int h, int ciclica) {
  printf("Resolucao: %d x %d\n", w, h);
  printf("MPI: %d processos | Distribuicao: %s\n\n", p, ciclica ? "ciclica" : "blocos");
  fflush(stdout);
}

// Exibe os tempos e o balanceamento de carga
static void mostra_tempos(int p, int sem_saida, double total, double minimo, double maximo, double soma, double espera_max, double coleta_max, double reordena, double escrita) {
  printf("  Tempo de Mandelbrot (calculo + coleta): %.6f s\n", total);
  printf("  Carga (%d processos): Tmin=%.6f s, Tmax=%.6f s, Tmedio=%.6f s, FLB=%.6f\n", p, minimo, maximo, soma / p, maximo > 0 ? (maximo - minimo) / maximo : 0);
  printf("  Comunicacao: espera=%.6f s, Gatherv=%.6f s, reordenacao=%.6f s\n", espera_max, coleta_max, reordena);
  if (!sem_saida)
    printf("  Tempo de escrita (BIN + PPM): %.6f s\n", escrita);
  putchar('\n');
  fflush(stdout);
}

// Salva a matriz BIN e a imagem PPM
static void salva_resultado(const int32_t *imagem, int w, int h, int limite, int ciclica, int cena, const char *nome) {
  char caminho[256];
  snprintf(caminho, sizeof caminho, OUTPUT_DIR "/%s_%d_%s.bin", ciclica ? "ciclica" : "blocos", cena + 1, nome);
  FILE *f = fopen(caminho, "wb");
  if (!f)
    falha("Erro ao criar BIN.");
  size_t gravados = fwrite(imagem, sizeof(int32_t), (size_t)w * h, f);
  int erro = fclose(f);
  if (gravados != (size_t)w * h || erro)
    falha("Erro ao gravar BIN.");
  snprintf(caminho, sizeof caminho, OUTPUT_DIR "/%s_%d_%s.ppm", ciclica ? "ciclica" : "blocos", cena + 1, nome);
  salva_ppm(imagem, w, h, limite, cena == 2, caminho);
}

// Lê e valida os argumentos do terminal
static void le_argumentos(int argc, char **argv, int config[5]) {
  if (argc < 2 || argc > 6 || (strcmp(argv[1], "blocos") && strcmp(argv[1], "ciclica")))
    falha("Uso: mandelbrot_mpi blocos|ciclica [largura altura [sem-saida [verificar]]]");
  config[2] = !strcmp(argv[1], "ciclica");
  if (argc == 3)
    falha("Informe largura e altura juntas.");
  if (argc >= 4) {
    config[0] = inteiro(argv[2]);
    config[1] = inteiro(argv[3]);
  }
  if (argc >= 5) {
    if (strcmp(argv[4], "sem-saida"))
      falha("Opcao esperada: sem-saida.");
    config[3] = 1;
  }
  if (argc == 6) {
    if (strcmp(argv[5], "verificar"))
      falha("Opcao esperada: verificar.");
    config[4] = 1;
  }
}

// Coordena o cálculo distribuído e a coleta dos resultados
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  // Largura, altura, estratégia, sem saída e verificação
  int config[5] = {4096, 4096, 0, 0, 0};
  if (!rank)
    le_argumentos(argc, argv, config);

  // Compartilha os parâmetros lidos pelo rank 0
  MPI_Bcast(config, 5, MPI_INT, 0, MPI_COMM_WORLD);

  int w = config[0], h = config[1], ciclica = config[2];
  if ((size_t)w * h > INT_MAX)
    falha("Gatherv limitado a INT_MAX pixels nesta implementacao.");
  int *contagens = aloca(p, sizeof(int)), *deslocamentos = aloca(p, sizeof(int));

  // Contagens e deslocamentos em pixels, distribui também as linhas restantes
  int deslocamento = 0;
  for (int r = 0; r < p; r++) {
    contagens[r] = (h / p + (r < h % p)) * w;
    deslocamentos[r] = deslocamento;
    deslocamento += contagens[r];
  }

  int nlinhas = contagens[rank] / w;

  // Blocos: linhas consecutivas
  // Cíclica: uma linha a cada p
  int primeira = ciclica ? rank : deslocamentos[rank] / w;
  int passo = ciclica ? p : 1;
  int32_t *local = aloca(contagens[rank], sizeof(int32_t));
  int32_t *reunido = !rank ? aloca((size_t)w * h, sizeof(int32_t)) : NULL;
  int32_t *imagem = !rank && ciclica ? aloca((size_t)w * h, sizeof(int32_t)) : reunido;
  if (!rank && !config[3])
    cria_pasta_saida();

  double regioes[4][4] = {
      {-2.0, 1.0, -1.5, 1.5},
      {-0.743643887 - 0.0015, -0.743643887 + 0.0015, 0.131825904 - 0.0015, 0.131825904 + 0.0015},
      {-2.0, 1.0, -1.5, 1.5},
      {-0.16070135 - 0.004, -0.16070135 + 0.004, 1.0375665 - 0.004, 1.0375665 + 0.004}};
  const int limites[] = {1000, 5000, 1000, 3000};
  const char *nomes[] = {"vista_completa", "vale_cavalos_marinhos", "mapa_de_custo", "zoom_espiral"};

  // Mantém a proporção do zoom em imagens retangulares
  double meia_altura = 0.0015 * (double)h / (double)w;
  regioes[1][2] = 0.131825904 - meia_altura;
  regioes[1][3] = 0.131825904 + meia_altura;
  const char *titulos[] = {"Vista completa", "Vale dos cavalos-marinhos", "Mapa de custo", "Zoom espiral"};
  if (!rank)
    mostra_configuracao(p, w, h, ciclica);

  for (int cena = 0; cena < 4; cena++) {
    if (!rank) {
      printf("Gerando %s...\n", titulos[cena]);
      fflush(stdout);
    }
    const double *a = regioes[cena];
    
    // Sincroniza antes de iniciar a medição
    MPI_Barrier(MPI_COMM_WORLD);
    double inicio = MPI_Wtime();
    for (int j = 0; j < nlinhas; j++)
      linha(local + (size_t)j * w, primeira + j * passo, w, h, limites[cena], a[0], a[1], a[2], a[3]);

    // Mede o cálculo antes da espera pelos outros processos
    double calculo = MPI_Wtime() - inicio;
    double t = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    double espera = MPI_Wtime() - t;
    t = MPI_Wtime();

    // Reúne as parcelas de tamanhos diferentes no rank 0
    MPI_Gatherv(local, contagens[rank], MPI_INT32_T, reunido, contagens, deslocamentos, MPI_INT32_T, 0, MPI_COMM_WORLD);
    double coleta = MPI_Wtime() - t;
    double reordena = 0;

    // Recoloca as linhas cíclicas na ordem da imagem
    if (!rank && ciclica) {
      t = MPI_Wtime();
      for (int r = 0; r < p; r++)
        for (int j = 0; j < contagens[r] / w; j++)
          memcpy(imagem + (size_t)(r + j * p) * w, reunido + deslocamentos[r] + (size_t)j * w, (size_t)w * sizeof(int32_t));
      reordena = MPI_Wtime() - t;
    }

    // Tempo total inclui espera, coleta e reordenação; exclui escrita e validação
    double total_local = MPI_Wtime() - inicio, total;
    MPI_Reduce(&total_local, &total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    double minimo, maximo, soma, espera_max, coleta_max;
    MPI_Reduce(&calculo, &minimo, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&calculo, &maximo, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&calculo, &soma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&espera, &espera_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&coleta, &coleta_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
      // Comparação exata das contagens com o sequencial
      if (config[4]) {
        int32_t *ref = mandelbrot(a[0], a[1], a[2], a[3], w, h, limites[cena]);
        if (memcmp(ref, imagem, (size_t)w * h * sizeof(int32_t)))
          falha("ERRO: matriz difere do sequencial.");
        free(ref);
        puts("  Validacao: OK.");
      }
      double escrita = 0;
      if (!config[3]) {
        t = MPI_Wtime();
        salva_resultado(imagem, w, h, limites[cena], ciclica, cena, nomes[cena]);
        escrita = MPI_Wtime() - t;
      }
      mostra_tempos(p, config[3], total, minimo, maximo, soma, espera_max, coleta_max, reordena, escrita);
    }
  }
  if (!rank && ciclica)
    free(imagem);
  free(reunido);
  free(local);
  free(contagens);
  free(deslocamentos);
  if (!rank)
    printf("Concluido.%s\n", config[3] ? "" : " Imagens salvas em: " OUTPUT_DIR);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
