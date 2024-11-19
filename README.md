Esse código MATLAB tem por objetivo simular um sistema de comunicação digital Massive MIMO que utiliza modulação QAM e multiplexação por divisão de frequência ortogonal (OFDM).

### Definição de parâmetros principais

- <b>SNR:</b> Define um vetor de valores para a Relação Sinal-Ruído (SNR) em decibéis, variando de -10 a 20 dB, para simular diferentes condições de ruído no canal de comunicação.
N_SNR: Calcula o comprimento do vetor SNR, ou seja, o número de valores de SNR que serão testados na simulação.

- <b>B:</b> Número de bits transmitidos por símbolo na modulação, definido aqui como 8 bits.
M_QAM: Calcula o tamanho da constelação para a modulação QAM (Quadrature Amplitude Modulation), que é $2^B$, resultando em uma constelação de 256 pontos (ou seja, 256-QAM).

- <b>M:</b> Número de subportadoras no sistema OFDM, ou seja, as portadoras independentes que serão usadas na multiplexação.
L: Comprimento do canal em número de amostras. Esse parâmetro representa a extensão do canal e a quantidade de amostras que ele retém do sinal.
K: Comprimento do prefixo cíclico em amostras. Para evitar interferência entre símbolos, o prefixo cíclico é utilizado e seu comprimento é igual ao comprimento do canal (L).
N: Comprimento total do bloco de dados, que é a soma das subportadoras (M) e do prefixo cíclico (K).

- <b>N_BLK:</b> Número de blocos de dados a serem transmitidos na simulação. Isso representa a quantidade de vezes que os dados serão enviados pelo canal durante a simulação.

- <b>N_EQ:</b> Define o número de equalizadores utilizados na simulação. A equalização é usada para compensar distorções no canal, e diferentes equalizadores podem ser testados para verificar seu desempenho.

- <b>N_FFT:</b> Número de pontos usados na Transformada Rápida de Fourier (FFT). A FFT é usada para converter o sinal do domínio do tempo para o domínio da frequência, essencial para o processamento em sistemas OFDM.

- <b>freq_N:</b> Frequências normalizadas associadas aos N_FFT pontos para a FFT, variando de -π a π.

- <b>freq_M:</b> Frequências associadas às subportadoras M, também variando de -π a π. Isso representa as frequências das subportadoras no domínio da frequência para o sistema OFDM.

Este conjunto de variáveis e parâmetros configura o sistema para simular a transmissão de dados em um canal de comunicação ruidoso, utilizando OFDM com modulação QAM, onde são analisados diferentes valores de SNR e utiliza-se FFT para o processamento no domínio da frequência.

### Separação de armazenamento

Esse trecho inicializa variáveis para simular a transmissão e recepção de dados em um sistema OFDM com diferentes tipos de prefixos.

<details>
    <summary><code>Conceito de prefixo</code></summary>

Aqui cabe uma pausa para explicação do conceito. "Prefixo", nesse contexto, refere-se a uma técnica usada para combater os efeitos da interferência entre símbolos (ISI, Inter-Symbol Interference) e a dispersão temporal causada por multipercursos no canal de comunicação. Existem diferentes tipos de prefixo, usaremos nesse trabalho os seguintes:

- <b>Prefixo Cíclico (Cyclic Prefix, CP):</b> Este é o tipo de prefixo mais comum em sistemas OFDM. Ele consiste em copiar uma parte final do símbolo OFDM e adicionar essa cópia no início do símbolo, criando um período de tempo extra antes do início do símbolo real. Facilita a equalização, pois transforma o canal em um sistema circular e evita a interferência entre símbolos ao criar uma zona de proteção contra o desvanecimento do canal.

- <b>Prefixo Zero (Zero Prefix, ZP):</b> Este prefixo adiciona uma sequência de zeros antes do símbolo OFDM em vez de copiar uma parte dele. Embora seja menos comum que o prefixo cíclico, o prefixo zero também pode reduzir a interferência entre símbolos ao fornecer uma janela onde os efeitos de interferência de múltiplos caminhos são reduzidos.
</details>

---

- <b>BER:</b> Matriz para armazenar a Taxa de Erro de Bits (Bit Error Rate, BER) para diferentes valores de SNR, número de equalizadores e tipos de prefixo (possivelmente prefixo cíclico, prefixo zero e outras variações). O tamanho da matriz é N_SNR (valores de SNR) x N_EQ (número de equalizadores) x 4 (tipos de prefixo).

- <b>A_cp:</b> Matriz de mapeamento para prefixo cíclico (cyclic prefix, CP). Essa matriz adiciona um prefixo cíclico ao bloco de dados OFDM.

- <b>zeros(K,M-K) eye(K):</b> Preenche as primeiras K linhas com zeros e adiciona uma matriz identidade de tamanho K.
    - <b>eye(M):</b> Matriz identidade de tamanho M, que representa as subportadoras principais.
    - <b>R_cp:</b> Matriz de remoção do prefixo cíclico. Após a transmissão, o prefixo cíclico é removido para recuperar os dados OFDM originais. A matriz seleciona as M subportadoras principais e ignora o prefixo cíclico.

<details>
  <summary><code>eye</code></summary>

Gera uma matriz identidade de tamanho especificado. A matriz identidade é uma matriz quadrada onde todos os elementos da diagonal principal são iguais a 1, enquanto todos os demais elementos são 0. Ela é comumente usada em operações de álgebra linear e processamento de sinais.

### Sintaxe
- `eye(n)`: Cria uma matriz identidade de tamanho `n x n`.
- `eye(m, n)`: Cria uma matriz de tamanho `m x n`, onde os elementos da diagonal principal são 1, e o restante é 0. Isso resulta em uma matriz identidade retangular se `m` e `n` forem diferentes.

### Exemplos
```matlab
eye(3)
% Resultado:
% 1 0 0
% 0 1 0
% 0 0 1

eye(2, 3)
% Resultado:
% 1 0 0
% 0 1 0
```

A matriz identidade é útil em operações onde se deseja manter os valores originais ao multiplicar vetores ou outras matrizes, já que, em álgebra linear, multiplicar por uma matriz identidade não altera o vetor ou matriz multiplicado.
</details>

---

- <b>A_zp:</b> Matriz de mapeamento para prefixo zero (zero prefix, ZP). Composta por uma matriz identidade eye(M) para as subportadoras e zeros(K,M) que adiciona K linhas de zeros no final, implementando o prefixo zero.

- <b>R_zp:</b> Matriz de remoção do prefixo zero. Remove o prefixo zero ao selecionar as M subportadoras principais do sinal recebido, ignorando os zeros.

- <b>x:</b> Matriz que armazena os dados transmitidos (possivelmente representados em diferentes prefixos) com tamanho N*N_BLK (total de dados transmitidos) x 4 (tipos de prefixo).

- <b>Px:</b> Vetor que armazena a potência dos dados transmitidos para cada tipo de prefixo. Essa potência será usada para normalizar o sinal e comparar o desempenho entre os prefixos.

- <b>y_p:</b> Matriz para armazenar os sinais recebidos após a transmissão pelo canal para cada bloco, tipo de prefixo, e valor de SNR. O tamanho é N x N_BLK x 4 x N_SNR.

- <b>m_mod_hat:</b> Matriz para armazenar os símbolos QAM demodulados (ou equalizados) após a recepção. Cada dimensão representa:

    - <b>M:</b> Número de subportadoras.
    - <b>N_BLK:</b> Número de blocos de dados.
    - <b>N_EQ:</b> Número de equalizadores.
    - <b>4:</b> Diferentes prefixos.
    - <b>N_SNR:</b> Diferentes valores de SNR.
    - <b>m_bit_hat:</b> Matriz para armazenar os bits estimados após a demodulação. Cada dimensão representa:

- <b>B\*M*N_BLK:</b> Número total de bits transmitidos.
    - <b>N_EQ:</b> Número de equalizadores.
    - <b>4:</b> Diferentes prefixos.
    - <b>N_SNR:</b> Diferentes valores de SNR.

Essas variáveis permitem simular a transmissão e recepção em diferentes condições de ruído, tipos de prefixo e configurações de equalização, registrando o desempenho (BER) para posterior análise.

---

### Modelagem do sistema de transmissão

Esse trecho está relacionado à geração e modulação de dados para transmissão em um sistema de comunicação digital utilizando modulação QAM (Quadrature Amplitude Modulation).

- <b>m_bit:</b> Cria um vetor de bits aleatórios com B*M*N_BLK elementos.

- <b>randi([0 1],B\*M\*N_BLK,1):</b> A função randi([0 1],n,1) gera n números aleatórios, que neste caso são bits (0 ou 1). O número total de bits gerados é B\*M\*N_BLK, que é calculado multiplicando o número de bits por símbolo (B), o número de subportadoras (M), e o número de blocos de dados (N_BLK).

- <b>m_mod:</b> Realiza a modulação QAM nos bits gerados.
A função qammod é usada para realizar a modulação QAM, que converte os bits de entrada em símbolos QAM.
O parâmetro 'InputType','bit' indica que a entrada para a modulação são bits (em vez de números inteiros, por exemplo).
M_QAM é o número de símbolos na constelação QAM, que foi definido anteriormente como 256 (ou seja, M_QAM = 2^B).

- <b>m_mod:</b> A função reshape é usada para reorganizar o vetor de símbolos modulados (m_mod) em uma matriz de tamanho M x N_BLK. M é o número de subportadoras, e N_BLK é o número de blocos de dados. Essa reorganização é importante para ajustar os dados modulados à estrutura do sistema OFDM, onde cada coluna da matriz representa um bloco de dados a ser transmitido, e cada linha representa uma subportadora.

<b>Detalhe:</b> Aqui, a variável m_mod está sendo modificada, não declarada novamente.

% Geração do ruído AWGN 
v = randn(N_BLK, 1) + 1i*randn(N_BLK, 1); 
% Ruído AWGN 
Pv = norm(v)^2 / length(v); 
% Potência média do ruído % Captura da saída com ruído adicionado 
y_ZF = H.' * x_ZF_normalized + v.'; 
% Saída ZF com ruído 
y_MF = H.' * x_MF_normalized + v.'; 
% Saída MF com ruído
% Exibir as saídas normalizadas com ruído 
disp('y_ZF (com ruído):'); 
disp(y_ZF); 
disp('y_MF (com ruído):'); 
disp(y_MF);