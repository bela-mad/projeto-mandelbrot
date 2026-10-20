/*
 * Geração paralela do conjunto de Mandelbrot com OpenMP.
 *
 * Comandos para compilar e executar (Windows):
 * gcc mandelbrot_openmp.c -o mandelbrot_openmp.exe -O3 -msse2 -mfpmath=sse -fopenmp -lm
 *
 * Exemplo de execução no PowerShell:
 * $env:OMP_NUM_THREADS = "4"
 * $env:OMP_DYNAMIC = "FALSE"
 * $env:OMP_SCHEDULE = "dynamic,1"
 * .\mandelbrot_openmp.exe
 *
 * O número de threads e o escalonamento podem ser alterados pelas variáveis OMP_NUM_THREADS e OMP_SCHEDULE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#ifdef MEDIR_CARGA
#include <float.h>
#endif

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#endif

// Diretório das imagens geradas
#define OUTPUT_DIR "saida"

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;

// Retorna o tempo monotônico em segundos
double tempo_atual() {
#ifdef _WIN32
    LARGE_INTEGER frequencia;
    LARGE_INTEGER contador;

    QueryPerformanceFrequency(&frequencia);
    QueryPerformanceCounter(&contador);

    return (double)contador.QuadPart / (double)frequencia.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
#endif
}

// Cria o diretório de saída, caso ainda não exista
void cria_pasta_saida() {
#ifdef _WIN32
    _mkdir(OUTPUT_DIR);
#else
    mkdir(OUTPUT_DIR, 0755);
#endif
}

// Calcula o número de iterações para cada ponto do plano complexo
int32_t *mandelbrot(double re_min, double re_max, double im_min, double im_max, int width, int height, int max_iter) {
    size_t total = (size_t)width * (size_t)height;
    int32_t *count = malloc(total * sizeof(int32_t));

    if (count == NULL) {
        fprintf(stderr, "Erro: nao foi possivel alocar memoria.\n");
        exit(EXIT_FAILURE);
    }

    // Distribui as linhas da imagem entre as threads (as variáveis declaradas no laço são privadas)
    // Cada thread escreve em posições distintas da matriz count
#ifdef MEDIR_CARGA
    double tempo_minimo = DBL_MAX;
    double tempo_maximo = 0.0;
    double soma_tempos = 0.0;
    int num_threads = 0;

    #pragma omp parallel default(none) \
        shared(count, re_min, re_max, im_min, im_max, width, height, max_iter) \
        reduction(min:tempo_minimo) reduction(max:tempo_maximo) \
        reduction(+:soma_tempos, num_threads)
    {
        // Sincroniza as threads antes de iniciar a medição
        #pragma omp barrier
        double inicio_thread = omp_get_wtime();

        // Permite medir a conclusão antes da espera pelas outras threads
        #pragma omp for schedule(runtime) nowait
#else
    #pragma omp parallel for default(none) \
        shared(count, re_min, re_max, im_min, im_max, width, height, max_iter) \
        schedule(runtime)
#endif
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

#ifdef MEDIR_CARGA
        double tempo_thread = omp_get_wtime() - inicio_thread;
        tempo_minimo = tempo_thread;
        tempo_maximo = tempo_thread;
        soma_tempos += tempo_thread;
        num_threads += 1;
    }

    // A barreira ao final da região parallel conclui a matriz e as reduções
    double fator = tempo_maximo > 0.0 ? (tempo_maximo - tempo_minimo) / tempo_maximo : 0.0;

    printf("  Carga (%d threads): Tmin=%.6f s, Tmax=%.6f s, Tmedio=%.6f s, FLB=%.6f\n",
           num_threads, tempo_minimo, tempo_maximo, soma_tempos / num_threads, fator);
#endif

    return count;
}

// Mantém o valor dentro dos limites mínimo e máximo
double limitar(double valor, double minimo, double maximo) {
    if (valor < minimo)
        return minimo;

    if (valor > maximo)
        return maximo;

    return valor;
}

// Realiza a interpolação linear entre duas cores
RGB interpola_cor(RGB a, RGB b, double t) {
    RGB resultado;

    resultado.r = (unsigned char)(a.r + (b.r - a.r) * t);
    resultado.g = (unsigned char)(a.g + (b.g - a.g) * t);
    resultado.b = (unsigned char)(a.b + (b.b - a.b) * t);

    return resultado;
}

// Aplica a paleta de cores Malha & Fluxo
RGB cor_malha_fluxo(double valor) {
    static const RGB paleta[] = {
        {11, 29, 58},   // #0b1d3a
        {18, 56, 110},  // #12386e
        {31, 111, 180}, // #1f6fb4
        {34, 182, 200}, // #22b6c8
        {240, 180, 41}, // #f0b429
        {255, 243, 209} // #fff3d1
    };

    const int num_cores = sizeof(paleta) / sizeof(paleta[0]);
    valor = limitar(valor, 0.0, 1.0);

    double posicao = valor * (num_cores - 1);
    int indice = (int)posicao;

    if (indice >= num_cores - 1)
        return paleta[num_cores - 1];

    double t = posicao - indice;

    return interpola_cor(paleta[indice], paleta[indice + 1], t);
}

// Aplica a paleta Inferno ao mapa de custo computacional
RGB cor_inferno(double valor) {
    static const RGB paleta[] = {
        {0, 0, 4},
        {40, 11, 84},
        {101, 21, 110},
        {159, 42, 99},
        {212, 72, 66},
        {245, 125, 21},
        {250, 193, 39},
        {252, 255, 164}
    };

    const int num_cores = sizeof(paleta) / sizeof(paleta[0]);
    valor = limitar(valor, 0.0, 1.0);

    double posicao = valor * (num_cores - 1);
    int indice = (int)posicao;

    if (indice >= num_cores - 1)
        return paleta[num_cores - 1];

    double t = posicao - indice;

    return interpola_cor(paleta[indice], paleta[indice + 1], t);
}

// Grava a matriz de contagens em arquivo binário, sem inverter a ordem das linhas
void salva_contagens(const int32_t *count, int width, int height, const char *path) {
    FILE *arquivo = fopen(path, "wb");

    if (arquivo == NULL) {
        fprintf(stderr, "Erro ao criar arquivo %s\n", path);
        exit(EXIT_FAILURE);
    }

    size_t total = (size_t)width * (size_t)height;
    size_t gravados = fwrite(count, sizeof(int32_t), total, arquivo);
    int erro = fclose(arquivo);

    if (gravados != total || erro != 0) {
        fprintf(stderr, "Erro ao gravar arquivo %s\n!", path);
        exit(EXIT_FAILURE);
    }

    printf("  Salvo: %s\n", path);
}

// Grava a representação colorida no formato PPM binário
void salva_colorida(const int32_t *count, int width, int height, int max_iter, const char *path) {
    FILE *arquivo = fopen(path, "wb");

    if (arquivo == NULL) {
        fprintf(stderr, "Erro ao criar arquivo %s\n!", path);
        return;
    }

    fprintf(arquivo, "P6\n%d %d\n255\n", width, height);

    double max_log = log1p((double)max_iter);

    for (int py = height - 1; py >= 0; py--) {
        for (int px = 0; px < width; px++) {
            int32_t iter = count[(size_t)py * width + px];
            RGB cor;

            if (iter >= max_iter) {
                cor.r = (unsigned char)(0.043 * 255.0);
                cor.g = (unsigned char)(0.114 * 255.0);
                cor.b = (unsigned char)(0.227 * 255.0);
            } else {
                double valor = log1p((double)iter) / max_log;
                cor = cor_malha_fluxo(valor);
            }

            fwrite(&cor, sizeof(RGB), 1, arquivo);
        }
    }

    fclose(arquivo);
    printf("  Salvo: %s\n", path);
}

// Grava o mapa de custo computacional no formato PPM binário
void salva_custo(const int32_t *count, int width, int height, int max_iter, const char *path) {
    FILE *arquivo = fopen(path, "wb");

    if (arquivo == NULL) {
        fprintf(stderr, "Erro ao criar arquivo %s\n!", path);
        return;
    }

    fprintf(arquivo, "P6\n%d %d\n255\n", width, height);

    for (int py = height - 1; py >= 0; py--) {
        for (int px = 0; px < width; px++) {
            int32_t iter = count[(size_t)py * width + px];
            double valor = (double)iter / (double)max_iter;
            RGB cor = cor_inferno(valor);

            fwrite(&cor, sizeof(RGB), 1, arquivo);
        }
    }
    fclose(arquivo);
    printf("  Salvo: %s\n", path);
}

int main(int argc, char *argv[]) {
    int width = 4096;
    int height = 4096;

    if (argc >= 2)
        width = atoi(argv[1]);

    if (argc >= 3)
        height = atoi(argv[2]);

    if (width < 2 || height < 2) {
        fprintf(stderr, "A resolução deve ser maior que 1 x 1.\n");
        return EXIT_FAILURE;
    }

    printf("Resolução: %d x %d\n\n", width, height);

    cria_pasta_saida();

    // Vista completa do conjunto
    printf("Gerando vista completa...\n");

    double inicio = tempo_atual();
    int32_t *full = mandelbrot(-2.0, 1.0, -1.5, 1.5, width, height, 1000);
    double fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Vista Completa): %.3f s\n", fim - inicio);

    inicio = tempo_atual();
    salva_contagens(full, width, height, OUTPUT_DIR "/mandelbrot_1_vista_completa.bin");
    salva_colorida(full, width, height, 1000, OUTPUT_DIR "/mandelbrot_1_vista_completa.ppm");
    fim = tempo_atual();

    printf("  Tempo de escrita (BIN + PPM): %.3f s\n", fim - inicio);
    free(full);

    // Ampliação da região conhecida como Vale dos Cavalos-Marinhos
    printf("Gerando vale dos cavalos-marinhos...\n");

    double cx = -0.743643887;
    double cy = 0.131825904;
    double half = 3.0e-3 / 2.0;

    inicio = tempo_atual();
    int32_t *sea = mandelbrot(cx - half, cx + half, cy - half, cy + half, width, height, 5000);
    fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Vale dos Cavalos-marinhos): %.3f s\n", fim - inicio);

    inicio = tempo_atual();
    salva_contagens(sea, width, height, OUTPUT_DIR "/mandelbrot_2_vale_cavalos_marinhos.bin");
    salva_colorida(sea, width, height, 5000, OUTPUT_DIR "/mandelbrot_2_vale_cavalos_marinhos.ppm");
    fim = tempo_atual();

    printf("  Tempo de escrita (BIN + PPM): %.3f s\n", fim - inicio);
    free(sea);

    // Mapa da distribuição do custo computacional
    printf("Gerando mapa de custo...\n");

    inicio = tempo_atual();
    int32_t *custo = mandelbrot(-2.0, 1.0, -1.5, 1.5, width, height, 1000);
    fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Mapa de Custo): %.3f s\n", fim - inicio);

    inicio = tempo_atual();
    salva_contagens(custo, width, height, OUTPUT_DIR "/mandelbrot_3_mapa_de_custo.bin");
    salva_custo(custo, width, height, 1000, OUTPUT_DIR "/mandelbrot_3_mapa_de_custo.ppm");
    fim = tempo_atual();

    printf("  Tempo de escrita (BIN + PPM): %.3f s\n", fim - inicio);
    free(custo);

    // Ampliação de uma região com formação espiral
    printf("Gerando zoom espiral...\n");

    double sx = -0.16070135;
    double sy = 1.0375665;
    double h2 = 0.004;

    inicio = tempo_atual();
    int32_t *spiral = mandelbrot(sx - h2, sx + h2, sy - h2, sy + h2, width, height, 3000);
    fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Zoom Espiral): %.3f s\n", fim - inicio);

    inicio = tempo_atual();
    salva_contagens(spiral, width, height, OUTPUT_DIR "/mandelbrot_4_zoom_espiral.bin");
    salva_colorida(spiral, width, height, 3000, OUTPUT_DIR "/mandelbrot_4_zoom_espiral.ppm");
    fim = tempo_atual();

    printf("  Tempo de escrita (BIN + PPM): %.3f s\n", fim - inicio);
    free(spiral);

    printf("\nConcluído. Imagens salvas na pasta: %s\n", OUTPUT_DIR);

    return 0;
}