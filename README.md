# Benchmark de Geração do Conjunto de Mandelbrot

Projeto desenvolvido para a disciplina **DEC107 – Processamento Paralelo**, do curso de Ciência da Computação da Universidade Estadual de Santa Cruz (UESC).

O projeto compara implementações do conjunto de Mandelbrot em C, analisando tempo de execução, speedup, eficiência e balanceamento de carga.

- **Etapa 1 — OpenMP:** paralelismo em memória compartilhada.
- **Etapa 2 — MPI:** paralelismo em memória distribuída, com distribuição por blocos e cíclica de linhas.

**Autoras:** Alícia Oliveira Araújo e Isabela Madureira Argolo

## Implementações

O projeto possui três versões principais:

- [`Etapa 1/mandelbrot_sequencial.c`](Etapa%201/mandelbrot_sequencial.c): implementação sequencial em C;
- [`Etapa 1/mandelbrot_openmp.c`](Etapa%201/mandelbrot_openmp.c): implementação paralela utilizando OpenMP;
- [`Etapa 2/mandelbrot_mpi.c`](Etapa%202/mandelbrot_mpi.c): implementação paralela utilizando MPI.

Os programas geram quatro representações do conjunto:

1. Vista completa;
2. Vale dos Cavalos-Marinhos;
3. Mapa de custo computacional;
4. Zoom em uma região com formação em espiral.

Os arquivos gerados ficam em `saida/` nas versões da etapa 1 e em `saida_mpi/` na versão MPI, dentro da pasta de execução. A implementação MPI mantém o cálculo por pixel, sem otimização por simetria.

## Requisitos

Os comandos abaixo são para o terminal PowerShell, que pode ser utilizado dentro do VS Code.

- GCC de 64 bits, com suporte a OpenMP para a etapa 1 e para os testes comparativos;
- Microsoft MPI Runtime e Microsoft MPI SDK para a etapa 2;
- comandos `gcc` e `mpiexec` disponíveis no terminal.

## Etapa 1 — Sequencial e OpenMP

Os comandos desta seção devem ser executados na pasta `Etapa 1`. A partir da raiz do projeto:

```powershell
cd ".\Etapa 1"
```

### Versão sequencial

No Windows com GCC:

```bash
gcc mandelbrot_sequencial.c -o mandelbrot_sequencial.exe -O3 -std=c11 -msse2 -mfpmath=sse -lm
```

Execução:

```powershell
.\mandelbrot_sequencial.exe
```

### Versão OpenMP

```bash
gcc mandelbrot_openmp.c -o mandelbrot_openmp.exe -O3 -std=c11 -msse2 -mfpmath=sse -fopenmp -lm
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

### Resolução

A versão OpenMP utiliza `4096 × 4096` por padrão, mas também permite informar a resolução pela linha de comando:

```powershell
.\mandelbrot_openmp.exe 2048 2048
```

Esse recurso foi utilizado principalmente nos testes de escalabilidade fraca.

### Balanceamento de carga

Para compilar a versão utilizada na medição do balanceamento entre as threads:

```bash
gcc mandelbrot_openmp.c -o mandelbrot_openmp_carga.exe -O3 -msse2 -mfpmath=sse -fopenmp -DMEDIR_CARGA -lm
```

A execução informa:

- menor tempo entre as threads (`Tmin`);
- maior tempo (`Tmax`);
- tempo médio (`Tmedio`);
- fator de balanceamento de carga (`FLB`).

### Testes da etapa 1

Os experimentos incluíram:

- execução sequencial;
- escalabilidade forte;
- escalabilidade fraca;
- comparação entre as políticas `static`, `dynamic` e `guided`;
- diferentes tamanhos de `chunk`;
- análise de balanceamento de carga.

Cada configuração foi executada três vezes, utilizando a mediana dos tempos para a análise.

### Saídas

A versão OpenMP gera:

- imagens no formato `.ppm`;
- matrizes de contagem no formato `.bin`.

A escrita dos arquivos é medida separadamente do tempo utilizado para o cálculo do Mandelbrot.

## Etapa 2 — MPI

Nesta etapa, as linhas da imagem são divididas entre processos. Foram implementadas duas estratégias:

- **Blocos:** cada processo calcula uma faixa de linhas contíguas.
- **Cíclica:** as linhas são distribuídas alternadamente entre os processos.

Cada processo calcula seus pixels, e o processo 0 reúne os resultados com `MPI_Gatherv`. Na distribuição cíclica, ele também reorganiza as linhas na ordem da imagem.

### Compilação

A partir da **raiz do repositório**, entre na pasta da etapa 2:

```powershell
cd ".\Etapa 2"
```

Compile com os caminhos padrão do Microsoft MPI SDK:

```powershell
gcc mandelbrot_mpi.c -O3 -std=c11 -msse2 -mfpmath=sse -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include" -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi -lm -o mandelbrot_mpi.exe
```

### Execução

Execute as estratégias com quatro processos:

```powershell
mpiexec -n 4 .\mandelbrot_mpi.exe blocos
mpiexec -n 4 .\mandelbrot_mpi.exe ciclica
```

O número após `-n` define a quantidade de processos.

### Resolução

A versão MPI utiliza `4096 × 4096` por padrão. Para informar outra resolução:

```powershell
mpiexec -n 4 .\mandelbrot_mpi.exe ciclica 2048 2048
```

### Balanceamento de carga e comunicação

A execução informa:

- tempo de cálculo com espera, coleta e reordenação;
- menor tempo entre os processos (`Tmin`), maior tempo (`Tmax`) e tempo médio (`Tmedio`);
- fator de balanceamento de carga, identificado como `FLB` no programa;
- tempos de espera, `Gatherv` e reordenação;
- tempo de escrita dos arquivos, quando habilitada.

A espera na barreira também reflete o desbalanceamento entre processos.

### Testes da etapa 2

Os experimentos incluíram:

- comparação entre as distribuições por blocos e cíclica;
- comparação com as versões sequencial e OpenMP;
- escalabilidade forte e fraca;
- análise de balanceamento de carga e custo de coleta;
- verificação de corretude das matrizes de contagem.

Os testes utilizaram a vista completa e o Vale dos Cavalos-Marinhos, com 1, 2, 4 e 8 processos ou threads. Cada configuração de desempenho foi executada três vezes após aquecimento, utilizando a mediana dos tempos para a análise.

### Saídas

A versão MPI gera:

- imagens no formato `.ppm`;
- matrizes de contagem no formato `.bin`.

Os arquivos ficam em `saida_mpi/`. A escrita é medida separadamente do tempo de cálculo com coleta.

## Nota de Transparência sobre o Uso de IA

Este projeto contou com o auxílio do **ChatGPT** como ferramenta de apoio na revisão do código, na preparação e correção do script de testes, na geração de gráficos e na revisão da escrita do relatório.

Todo o conteúdo produzido ou modificado com auxílio da ferramenta foi revisado, testado e validado pelas autoras. As autoras assumem total responsabilidade pela correção lógica do código, precisão dos resultados apresentados e integridade acadêmica do material entregue.
