# Benchmark de Geração do Conjunto de Mandelbrot

Projeto desenvolvido para a disciplina **DEC107 – Processamento Paralelo**, do curso de Ciência da Computação da Universidade Estadual de Santa Cruz (UESC).

O objetivo desta etapa é comparar uma implementação sequencial do conjunto de Mandelbrot com uma versão paralela utilizando **OpenMP**, analisando tempo de execução, speedup, eficiência, políticas de escalonamento e balanceamento de carga.

**Autoras:** Alícia Oliveira Araújo e Isabela Madureira Argolo

## Implementações

O projeto possui duas versões principais:

- `mandelbrot_sequencial.c`: implementação sequencial em C;
- `mandelbrot_openmp.c`: implementação paralela utilizando OpenMP.

As duas versões geram quatro representações do conjunto:

1. Vista completa;
2. Vale dos Cavalos-Marinhos;
3. Mapa de custo computacional;
4. Zoom em uma região com formação em espiral.

Os arquivos gerados são armazenados na pasta `saida/`.

## Compilação

### Versão sequencial

No Windows com GCC:

```bash
gcc mandelbrot_sequencial.c -o mandelbrot_sequencial.exe -O3 -lm
```

Execução:

```powershell
.\mandelbrot_sequencial.exe
```

### Versão OpenMP

```bash
gcc mandelbrot_openmp.c -o mandelbrot_openmp.exe -O3 -msse2 -mfpmath=sse -fopenmp -lm
```

Exemplo de execução no PowerShell:

```powershell
$env:OMP_NUM_THREADS = "4"
$env:OMP_DYNAMIC = "FALSE"
$env:OMP_SCHEDULE = "dynamic,1"

.\mandelbrot_openmp.exe
```

O número de threads pode ser alterado por `OMP_NUM_THREADS` e a política de escalonamento por `OMP_SCHEDULE`.

Exemplos:

```powershell
$env:OMP_SCHEDULE = "static,1"
$env:OMP_SCHEDULE = "dynamic,8"
$env:OMP_SCHEDULE = "guided,32"
```

## Resolução

A versão OpenMP utiliza `4096 × 4096` por padrão, mas também permite informar a resolução pela linha de comando:

```powershell
.\mandelbrot_openmp.exe 2048 2048
```

Esse recurso foi utilizado principalmente nos testes de escalabilidade fraca.

## Balanceamento de carga

Para compilar a versão utilizada na medição do balanceamento entre as threads:

```bash
gcc mandelbrot_openmp.c -o mandelbrot_openmp_carga.exe -O3 -msse2 -mfpmath=sse -fopenmp -DMEDIR_CARGA -lm
```

A execução informa:

- menor tempo entre as threads (`Tmin`);
- maior tempo (`Tmax`);
- tempo médio (`Tmedio`);
- fator de balanceamento de carga (`FLB`).

## Testes realizados

Os experimentos incluíram:

- execução sequencial;
- escalabilidade forte;
- escalabilidade fraca;
- comparação entre as políticas `static`, `dynamic` e `guided`;
- diferentes tamanhos de `chunk`;
- análise de balanceamento de carga.

Cada configuração foi executada três vezes, utilizando a mediana dos tempos para a análise.

## Saídas

A versão paralela gera:

- imagens no formato `.ppm`;
- matrizes de contagem no formato `.bin`.

A escrita dos arquivos é medida separadamente do tempo utilizado para o cálculo do Mandelbrot.

## Nota de Transparência sobre o Uso de IA

Este projeto contou com o auxílio do **ChatGPT** como ferramenta de apoio na revisão do código, estruturação dos testes, geração de gráficos e revisão da escrita do relatório.

Todo o conteúdo produzido ou modificado com auxílio da ferramenta foi revisado, testado e validado pelas autoras. As autoras assumem total responsabilidade pela correção lógica do código, precisão dos resultados apresentados e integridade acadêmica do material entregue.

**Alícia Oliveira Araújo e Isabela Madureira Argolo – 22/09/2026**
