close all
clear all
clc
pkg load control


% 1º: Geração da população inicial com valores aleatórios
function Populacao = InicializaPopulacao(tam_populacao,tam_cromossomos,limite_min,limite_max,limite_max_KiKd)
  Populacao = zeros(tam_populacao,tam_cromossomos); %inicializa a matriz de população aleatória: 16 linhas e 3 colunas
  Populacao(:,1) = rand(tam_populacao,1)*(limite_max - limite_min)+limite_min;  %gera valores aleatórios de até 100 para Kp
  Populacao(:,2) = rand(tam_populacao,1)*limite_max_KiKd; %gera valores aleatórios de até 500 para Ki
  Populacao(:,3) = rand(tam_populacao,1)*limite_max_KiKd; %gera valores aleatórios de até 500 para Kd

end

function IAE = Simula_Avalia(Kp,Ki,Kd) %Função desempenho a ser otimizada (planta do sistema)
  m = 0.007;
  R = 0.05;
  L = 0.85;
  d = 0.092;
  g = -9.81;
  J = 7*10^-6;
s = tf('s');
BallBeam = tf(-m*g*d/L/(m+(J/R^2))/s^2);
Controlador = pid(Kp,Ki,Kd);
MF = feedback(Controlador*BallBeam,1);
setpoint = 0.2;
t = 0:0.01:40;  #antes estava 0:0.01:180
y = step(setpoint*MF,t);

%para o cálculo do valor absoluto do erro: IAE

IAE=0;
  for i = 1:length(t)
    erro = setpoint - y(i);
    IAE = IAE + abs(erro);
  end

end


function Fitness = AvaliaPopulacao(Populacao) %avalia cada membro da população, suas aptidões
  tam_populacao = size(Populacao,1);
  Fitness = zeros(tam_populacao, 1); %ele cria a matriz das aptidões dos invidíduos
  disp("\nAptidões dos indivíduos:\n")
  for i = 1:tam_populacao;
    Kp = Populacao(i, 1);
    Ki = Populacao(i, 2);
    Kd = Populacao(i, 3);
    Fitness(i) = Simula_Avalia(Kp,Ki,Kd);
    fprintf('Indivíduo %d - Kp: %f, Ki: %f, Kd: %f, IAE: %f\n', i, Kp, Ki, Kd, Fitness(i));
  end
end



#FUNÇÃO DE SELEÇÃO DOS PAIS: MÉTODO DE TORNEIO RETORNANDO 2 PAIS
function [melhores_individuos,indices_melhores] = SelecionaMelhorIndividuo_Torneio(Populacao,Fitness,tamanho_torneio)
  #temos de percorrer a população e selecionar aleatoriamente os indivíduos (estou fazendo isso pelos índices deles),
  %depois comparar as suas aptidões
  %o que tiver menor aptidão será o pai

  tam_populacao = size(Populacao,1);
  melhores_individuos = zeros(2, size(Populacao, 2)); %cria uma matriz 2x2 que armazenará os genes Kp,Ki e Kd dos melhores indivíduos selecionados
  indices_melhores = zeros(2,1); %isso cria um vetor de duas linhas e uma coluna contendo os índices do melhores indivíduos selecionados

  for i = 1:2; #para pegar os dois melhores indivíduos que serão os pais (um de cada torneio)
  indicesDoTorneio = randi(tam_populacao,tamanho_torneio,1); %seleciona índices aleatórios da população para o torneio
  aptidoesCompetidores = Fitness(indicesDoTorneio); %calcula as aptidões dos competidores do torneio
  [~, ind_melhor] = min(aptidoesCompetidores); %pega a melhor aptidão dos competidores (o melhor é o que tiver menor IAE)
  indices_melhores(i) = indicesDoTorneio(ind_melhor); #retorna o indice do melhor individuo
  melhores_individuos(i,:) = Populacao(indices_melhores(i), :); #retorna o melhor indivíduo
  endfor
end



#FUNÇÃO DE SELEÇÃO DOS PAIS: MÉTODO DA ROLETA RETORNANDO 2 PAIS
function [melhores_individuos2, indices_melhores2] = SelecionaMelhorIndividuo_Roleta(Populacao,Fitness)
  #Na seleção por roleta, cada indivíduo da população é representado na roleta proporcionalmente ao seu valor de aptidão.
  #Assim, para indivíduos com alta aptidão é dada uma porção maior da roleta, enquanto aos indivíduos de aptidão mais baixa,
  #é dada uma porção relativamente menor.

  #A representação das aptidões na roleta é feita da seguinte forma:
  #IAE(i)/(somatório dos IAEs)

  #depois é criada a roleta com as probabilidades dos indivíduos da população, depois a roleta é girada (gerando um número aleatório)
  #e selecionando os indivíduos específicos aonde ela parou em cada giro)
  #Ela será girada duas vezes para pegar os dois pais.

  aptidaoTotal = sum(Fitness); #calcula a aptidão total (somatório dos IAEs da população)
  probabilidades_dos_individuos = (Fitness/aptidaoTotal); #calcula a probabilidade de cada individuo ser selecionado na roleta

  #criação da roleta onde as posições são as probabilidades acumuladas dos indivíduos até aquela posição,
  #dessa forma, a seleção dos indivíduos fica mais simples posteriormente
  roleta = cumsum(probabilidades_dos_individuos);

  numero_de_pais = 2;

  #inicialização das matrizes que armazenarão os resultados
  melhores_individuos2 = zeros(numero_de_pais, size(Populacao, 2));
  indices_melhores2 = zeros(numero_de_pais,1);

  #looping para gerar os números aleatórios de 0 a 1 (giro da roleta) e selecionar os pais
  for i = 1:numero_de_pais; #gira duas vezes a roleta para pegar um pai em cada giro
    while true
      x = rand(); #giro da roleta que gera um valor de aleatório de 0 a 1
      for j = 1:length(roleta);
        if x <= roleta(j);
          if ~ismember(j, indices_melhores2); #isso garante que o mesmo pai não será selecionado novamente
            melhores_individuos2(i,:) = Populacao(j,:); #seleciona os pais
            indices_melhores2(i) = j; #pega os índices do pais
          break;
          endif
        endif
      endfor
        if indices_melhores2(i) != 0;
          break;
        endif
    endwhile
  endfor
end


#FUNÇÃO DE CROSSOVER 1: CROSSOVER DE DOIS PONTOS
function [filho_C1,filho_C2] = CrossoverDosPais_DoisPontos(pai1,pai2)
  #passamos o tamanho dos genes dos pais e selecionamos aleatoriamente os pontos de corte nos dois (que possuem o mesmo tamanho)
  #inicializamos os filhos como tendo os genes dos pais, depois fazemos a troca entre os filhos dos genes cortados anteriormente
  #e teremos então os novos indivíduos, os filhos

  tam_cromossomos = length(pai1); #passando o tamanho do pai

  #verificando os pontos de corte
  PrimeiroPonto = randi([1 tam_cromossomos-1]); #gera valores aleatórios para o primeiro ponto de corte dentro do range de valores dos genes
  SegundoPonto = randi([PrimeiroPonto+1 tam_cromossomos]); #gera valores aleatórios para o segundo ponto de corte dentro do range de valores dos genes

  #para gerar os filhos, primeiro nós copiamos os valores dos pais
  filho_C1 = pai1;
  filho_C2 = pai2;

  #agora realizamos as trocas dos genes, gerando os filhos:
  filho_C1(PrimeiroPonto:SegundoPonto) = pai2(PrimeiroPonto:SegundoPonto); #primeiro filho recebe os genes do segundo pai
  filho_C2(PrimeiroPonto:SegundoPonto) = pai1(PrimeiroPonto:SegundoPonto); #segundo filho recebe os genes do primeiro pai

end


#SEGUNDO MÉTODO DE CROSSOVER: CROSSOVER DE 1 PONTO
function [filho_CB1,filho_CB2] = CrossoverDosPais_UmPonto(pai1,pai2)

  tam_cromossomos = length(pai1); #passamos o tamanho dos cromossomos dos pais, como ambos possuem o mesmo tamanho, só passo 1

  UnicoCorte = randi([1 tam_cromossomos-1]); #gera o ponto de corte nos genes dos pais

  #incializa os filhos com os genes dos pais
  filho_CB1 = pai1;
  filho_CB2 = pai2;

  #reliza a troca dos genes
  filho_CB1(UnicoCorte) = pai2(UnicoCorte);
  filho_CB2(UnicoCorte) = pai1(UnicoCorte);

end


#MÉTODO DE MUTAÇÃO: MUTAÇÃO UNIFORME SOMANDO 1 NO GENE DO INDIVÍDUO
function mutacao = MutacaoUniformeDosFilhos(individuo,taxaDeMutacao,limite_min, limite_max, limite_max_KiKd)

  mutacao = individuo; #passa o valor do individuo para inicializar a mutacao
  numeroDeGenes = length(individuo); #passa o tamanho dos genes do individuo

  for i = 1:numeroDeGenes; #iterando sobre o tamanho dos genes do indivíduo

    if rand() < taxaDeMutacao; #gera um número aleatório de 0 a 1 e compara com a taxa de mutação
      if i == 1 #analisa o gene Kp
        incremento = mutacao(i)*0.01; #valor 1 em porcentagem
        if mutacao(i) + incremento > limite_max;
          mutacao(i) = limite_max; #forço o gene ter o valor do limite max da população para o Kp
        else
          mutacao(i) = mutacao(i) + incremento; #aplico a mutação se estiver tudo certo
        endif
      else #para Ki e Kd que possuem o mesmo valor de limite
        incremento = mutacao(i)*0.01;
        if mutacao(i) + incremento > limite_max_KiKd;
          mutacao(i) = limite_max_KiKd; #forço o gene ter o valor do limite max da população para o KiKd
        else
          mutacao(i) = mutacao(i) + incremento; #aplico a mutação se estiver tudo certo
        endif
      endif
      fprintf("Gene %d mutado em +1%%, seu novo valor é: %f\n",i,mutacao(i));

    else  # rand() > taxaDeMutacao;

      if i == 1 %analisa o gene Kp
        decremento = mutacao(i)*0.01;
        if mutacao(i) - decremento < limite_min;
          mutacao(i) = limite_min; #força o gene Kp ter o valor do limite mínimo da população
        else
          mutacao(i) = mutacao(i) - decremento; #aplico o decremento no gene Kp (mutação)
        endif
      else  #para os genes Ki e Kd
        decremento = mutacao(i)*0.01;
        if mutacao(i) - decremento < limite_min;
          mutacao(i) = limite_min;
        else
          mutacao(i) = mutacao(i) - decremento;
        endif
      endif
      fprintf("Gene %d mutado em -1%%, seu novo valor é: %f\n",i,mutacao(i));
    endif

  endfor

end



#SELEÇÃO PRINCIPAL DO ALGORITMO GENÉTICO QUE UNE TODAS AS FUNÇÕES
function [AptidaoMaxima, AptidaoMinima, Media_das_Aptidoes, melhor_individuo_geral, melhor_aptidao_geral] = AlgoritmoGenetico(tam_populacao,tam_cromossomos,limite_min,limite_max,limite_max_KiKd,tamanho_torneio,taxaDeMutacao,NumeroDeGeracoes)

  disp("\n\tPopulação aleatória:\n")
  disp("     Kp:\tKi:\t Kd:")
  Populacao = InicializaPopulacao(tam_populacao,tam_cromossomos,limite_min,limite_max,limite_max_KiKd);
  disp(Populacao)

  Kp=1;Ki=1;Kd=1;
  fprintf("\n IAE para Kp = %d, Ki = %d, Kd = %d: %f\n", Kp, Ki, Kd, Simula_Avalia(Kp, Ki, Kd))

  Fitness = AvaliaPopulacao(Populacao); #avalia a populacao

  #inicializando as variáveis abaixo para pegar os melhores indivíduos e suas aptidões em cada geração:
  melhores_aptidoes = zeros(NumeroDeGeracoes,1); #pega as melhores aptidôes de cada geração
  melhores_individuos_geracao = zeros(NumeroDeGeracoes, tam_cromossomos); #pega os melhores indivíduos de cada geração (Kp, Ki e Kd)

  melhor_aptidao_geral = inf; #começa com infinito
  melhor_individuo_geral = zeros(1, tam_cromossomos); #guarda o melhor individuo das gerações

  for geracao = 1:NumeroDeGeracoes
    disp("\n")
    disp(['Geração ', num2str(geracao)]);

    #realizamos a seleção,crossover e mutação

    #seleção:
    [melhores_individuos, indices_melhores] = SelecionaMelhorIndividuo_Torneio(Populacao, Fitness, tamanho_torneio); #seleção torneio
    #[melhores_individuos2,indices_melhores2] = SelecionaMelhorIndividuo_Roleta(Populacao,Fitness); #seleção roleta

    #crossover dois pontos:  C1 e C2
    #[filho_C1,filho_C2] = CrossoverDosPais_DoisPontos(melhores_individuos(1,:),melhores_individuos(2,:)); # passando a seleção por torneio
    #[filho_C1,filho_C2] = CrossoverDosPais_DoisPontos(melhores_individuos2(1,:),melhores_individuos2(2,:)); # passando a seleção por roleta

    #crossover um ponto:    CB1 e CB2
    [filho_CB1,filho_CB2] = CrossoverDosPais_UmPonto(melhores_individuos(1,:),melhores_individuos(2,:)); #passando a seleção por torneio
    #[filho_CB1,filho_CB2] = CrossoverDosPais_UmPonto(melhores_individuos2(1,:),melhores_individuos2(2,:)); #passando a seleção por roleta


    #mutação:
    #mutacao_filho_C1 = MutacaoUniformeDosFilhos(filho_C1,taxaDeMutacao, limite_min, limite_max, limite_max_KiKd); #após o crossover de dois pontos
    #disp("\nMutação do primeiro filho: ")
    #disp(mutacao_filho_C1)

    #mutacao_filho_C2 = MutacaoUniformeDosFilhos(filho_C2,taxaDeMutacao, limite_min, limite_max, limite_max_KiKd); #após o crossover de dois pontos
    #disp("\nMutação do segundo filho: ")
    #disp(mutacao_filho_C2)

    mutacao_filho_CB1 = MutacaoUniformeDosFilhos(filho_CB1,taxaDeMutacao, limite_min, limite_max, limite_max_KiKd); #após o crossover de um ponto
    disp("\nMutação do primeiro filho: ")
    disp(mutacao_filho_CB1)

    mutacao_filho_CB2 = MutacaoUniformeDosFilhos(filho_CB2,taxaDeMutacao, limite_min, limite_max, limite_max_KiKd);
    disp("\nMutação do segundo filho: ")
    disp(mutacao_filho_CB2)

    # Obtenha o número de indivíduos a serem substituídos
    num_individuos_substituidos = 2;

    # Verifique o tamanho das mutações
    if size([mutacao_filho_CB1; mutacao_filho_CB2], 1) == num_individuos_substituidos
      # Ordene a população com base no Fitness e obtenha os índices dos piores indivíduos
      [~, indices_piores] = sort(Fitness, 'ascend');
      indices_piores = indices_piores(1:num_individuos_substituidos);

      # Substituição dos piores indivíduos pela nova geração
      Populacao(indices_piores, :) = [mutacao_filho_CB1; mutacao_filho_CB2];
    else
      error('Número de indivíduos para substituição e novos indivíduos não coincidem.');
    endif

    Fitness = AvaliaPopulacao(Populacao); #Avalia a nova populacao

    #Atualiza os melhores indivíduos e suas aptidôes a cada geração:
    [melhor_valor,indice_melhor] = min(Fitness);
    melhores_aptidoes(geracao) = melhor_valor;
    melhores_individuos_geracao(geracao, :) = Populacao(indice_melhor, :);

    #Atualiza a melhor aptidão geral entre todas as gerações:
    if melhor_valor < melhor_aptidao_geral
      melhor_aptidao_geral = melhor_valor;
      melhor_individuo_geral = Populacao(indice_melhor, :);
    endif

    disp("\nPopulação após geração:");
    disp("     Kp:\tKi:\t Kd:")
    disp(Populacao);
    disp("\nAptidões após geração:");
    disp(Fitness);

  endfor

  % Exibe o melhor indivíduo da última geração (será os nossos valores de Kp, Ki e Kd -> É O RESULTADO ESPERADO DO ALGORITMO GENÉTICO
  [melhor_valor, indice_melhor] = min(Fitness);
  disp("\n");
  disp(['MELHOR INDIVÍDUO DO ALGORITMO GENÉTICO APÓS A ÚLTIMA GERAÇÃO: Kp = ', num2str(Populacao(indice_melhor, 1)), ', Ki = ', num2str(Populacao(indice_melhor, 2)), ', Kd = ', num2str(Populacao(indice_melhor, 3)), ', IAE = ', num2str(melhor_valor)]);

  disp("\n");
  #exibindo o melhor indivíduo geral de todas as gerações:
  disp(['Melhor indivíduo de todas as gerações: Kp = ',num2str(melhor_individuo_geral(1)), ', Ki = ', num2str(melhor_individuo_geral(2)), ', Kd = ', num2str(melhor_individuo_geral(3)), ', IAE = ', num2str(melhor_aptidao_geral)]);

  disp("\n");
  disp('Vetor com as melhores aptidões de cada geração (IAE): ')
  disp(melhores_aptidoes)

  #cálculo do mínimo, máximo e média das melhores aptidôes de geração:
  AptidaoMaxima = max(melhores_aptidoes);
  AptidaoMinima = min(melhores_aptidoes);
  Media_das_Aptidoes = mean(melhores_aptidoes);

  disp("\n");
  disp(['Aptidão Máxima das melhores aptidões de cada geração: IAE = ', num2str(AptidaoMaxima)]);
  disp("\n");
  disp(['Aptidão Mínima das melhores aptidões de cada geração: IAE = ', num2str(AptidaoMinima)]);
  disp("\n");
  disp(['Média das melhores aptidões de cada geração: IAE = ', num2str(Media_das_Aptidoes)]);

end

tam_populacao = 20;  %tamanho da população
tam_cromossomos = 3;  %tamanho dos cromossomos (genes) Kp,Ki,Kd
limite_min = 0.5;
limite_max = 10; #estava 50
limite_max_KiKd = 3; #estava 10
tamanho_torneio = 4;
taxaDeMutacao = 0.1; #estava 0.5
NumeroDeGeracoes = 5; #estava 25
num_testes = 5;  % número de testes que queremos realizar

#AlgoritmoGenetico(tam_populacao,tam_cromossomos,limite_min,limite_max,limite_max_KiKd,tamanho_torneio,taxaDeMutacao,NumeroDeGeracoes)


#------------ DAQUI PARA BAIXO OS CÓDIGOS SÃO REFERENTES AO LOOPING DE TESTES EXAUSTIVOS FEITOS NO ALGORITMO -----------------

#inicializa os vetores para armazenar as aptidões máxima, mínima e a média de cada teste
#AptidoesMaximas_dosTestes = zeros(1, num_testes);
#AptidoesMinimas_dosTestes = zeros(1, num_testes);
#MediaDasAptidoes_dosTestes = zeros(1, num_testes);

#inicializa uma matriz que armazenará as aptidões máxima, mínima e a média de cada teste
resultado_dos_testes = zeros(num_testes, 3);

#para armazenar o melhor individúdo dentre todos os testes
melhor_individuo_dentre_TodosOstestes = zeros(1, tam_cromossomos);
melhor_aptidao_dentre_TodosOsTestes = inf; #começa infinito para que ter certeza de que qualquer aptidão será menor

disp("\n\n");
for teste = 1:num_testes
    disp(['----------------------------------------- Início do Teste ', num2str(teste), ' -----------------------------------------']);

    [AptidaoMaxima, AptidaoMinima, Media_das_Aptidoes, melhor_individuo_geral, melhor_aptidao_geral] = AlgoritmoGenetico(tam_populacao, tam_cromossomos, limite_min, limite_max, limite_max_KiKd, tamanho_torneio, taxaDeMutacao, NumeroDeGeracoes);

    #armazena a média, o máximo e o mínimo de cada teste de acordo com as variáveis das gerações
    #AptidoesMaximas_dosTestes(teste) = AptidaoMaxima;
    #AptidoesMinimas_dosTestes(teste) = AptidaoMinima;
    #MediaDasAptidoes_dosTestes(teste) = Media_das_Aptidoes;
    resultado_dos_testes(teste, :) = [AptidaoMaxima, AptidaoMinima, Media_das_Aptidoes];

    #verifica se o melhor indivíduo do teste em questão é melhor de todos
    if melhor_aptidao_geral < melhor_aptidao_dentre_TodosOsTestes
      melhor_aptidao_dentre_TodosOsTestes = melhor_aptidao_geral;
      melhor_individuo_dentre_TodosOstestes = melhor_individuo_geral;
    endif

    disp("\n");
    disp(['\n----------------------------------------- Fim do Teste ', num2str(teste), ' -----------------------------------------']);
    disp("\n");
end


#exibindo os resultados finais gerais de todos os testes
disp("\n\n\n");
disp("Resultados dos testes feitos, em forma de matriz: \n");
#disp(["Aptidão Máxima de cada Teste respectivamente: ", num2str(AptidoesMaximas_dosTestes)]);
#disp(["Aptidão Mínima de cada Teste respectivamente: ", num2str(AptidoesMinimas_dosTestes)]);
#disp(["Média das Aptidões de cada Teste respectivamente: ", num2str(MediaDasAptidoes_dosTestes)]);
disp("1ª Coluna - Aptidão Máxima de cada Teste respectivamente");
disp("2ª Coluna - Aptidão Mínima de cada Teste respectivamente");
disp("3ª Coluna - Média das Aptidões de cada Teste respectivamente");
disp("\n");
disp(resultado_dos_testes);
disp("\n\n");

#mostra o melhor indivíduo dentre todos os testes
disp(["MELHOR INDIVÍDUO DENTRE TODOS OS TESTES: Kp = ", num2str(melhor_individuo_dentre_TodosOstestes(1)), ", Ki = ", num2str(melhor_individuo_dentre_TodosOstestes(2)), ", Kd = ", num2str(melhor_individuo_dentre_TodosOstestes(3)), ", IAE = ", num2str(melhor_aptidao_dentre_TodosOsTestes)]);
disp("\n\n");


