/*
* Geração sequencial de imagens do conjunto de Mandelbrot
*
* Comandos para compilar e executar no Windows:
*
* gcc mandelbrot_sequencial.c -o mandelbrot_sequencial.exe -O3 -lm
* .\mandelbrot_sequencial.exe
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

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
        {252, 255, 164}};

    const int num_cores = sizeof(paleta) / sizeof(paleta[0]);
    valor = limitar(valor, 0.0, 1.0);
    double posicao = valor * (num_cores - 1);
    int indice = (int)posicao;

    if (indice >= num_cores - 1)
        return paleta[num_cores - 1];

    double t = posicao - indice;

    return interpola_cor(paleta[indice], paleta[indice + 1], t);
}

// Grava a representação colorida no formato PPM binário
void salva_colorida(const int32_t *count, int width, int height, int max_iter, const char *path) {
    FILE *arquivo = fopen(path, "wb");

    if (arquivo == NULL) {
        fprintf(stderr, "Erro ao criar arquivo %s\n", path);
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
        fprintf(stderr, "Erro ao criar arquivo %s\n", path);
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

int main() {
    cria_pasta_saida();

    // Vista completa do conjunto
    printf("Gerando vista completa...\n");

    double inicio = tempo_atual();
    int32_t *full = mandelbrot(-2.0, 1.0, -1.5, 1.5, 4096, 4096, 1000);
    double fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Vista Completa): %.3f s\n", fim - inicio);
    salva_colorida(full, 4096, 4096, 1000, OUTPUT_DIR "/mandelbrot_1_vista_completa.ppm");
    free(full);

    // Ampliação da região conhecida como Vale dos Cavalos-Marinhos
    printf("Gerando vale dos cavalos-marinhos...\n");

    double cx = -0.743643887;
    double cy = 0.131825904;
    double half = 3.0e-3 / 2.0;

    inicio = tempo_atual();
    int32_t *sea = mandelbrot(cx - half, cx + half, cy - half, cy + half, 4096, 4096, 5000);
    fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Vale dos Cavalos-marinhos): %.3f s\n", fim - inicio);
    salva_colorida(sea, 4096, 4096, 5000, OUTPUT_DIR "/mandelbrot_2_vale_cavalos_marinhos.ppm");
    free(sea);

    // Mapa da distribuição do custo computacional
    printf("Gerando mapa de custo...\n");

    inicio = tempo_atual();
    int32_t *custo = mandelbrot(-2.0, 1.0, -1.5, 1.5, 4096, 4096, 1000);
    fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Mapa de Custo): %.3f s\n", fim - inicio);
    salva_custo(custo, 4096, 4096, 1000, OUTPUT_DIR "/mandelbrot_3_mapa_de_custo.ppm");
    free(custo);

    // Ampliação de uma região com formação espiral
    printf("Gerando zoom espiral...\n");

    double sx = -0.16070135;
    double sy = 1.0375665;
    double h2 = 0.004;

    inicio = tempo_atual();
    int32_t *spiral = mandelbrot(sx - h2, sx + h2, sy - h2, sy + h2, 4096, 4096, 3000);
    fim = tempo_atual();

    printf("  Tempo de Mandelbrot (Zoom Espiral): %.3f s\n", fim - inicio);
    salva_colorida(spiral, 4096, 4096, 3000, OUTPUT_DIR "/mandelbrot_4_zoom_espiral.ppm");
    free(spiral);
    printf("\nConcluido. Imagens salvas na pasta: %s\n", OUTPUT_DIR);

    return 0;
}