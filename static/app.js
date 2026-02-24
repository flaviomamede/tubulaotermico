(function () {
  'use strict';

  const API_OTIMIZAR_CLOUD = '';
  const API_OTIMIZAR = API_OTIMIZAR_CLOUD || (window.location.origin + '/api/otimizar');
  const API_CURVA = (window.location.origin + '/api/curva');

  const NOMES_PARAMS = ['dT_adi1', 'dT_adi2', 'tau1', 'beta1', 'tau2', 'beta2', 'k_rel', 'alpha1', 'alpha2'];
  const DEFAULT_CHUTE = [45.0, 40.0, 10.0, 3.0, 25.0, 1.5, 2.9, 0.0040, 0.0030];
  const DEFAULT_BOUNDS_INF = [5.0, 5.0, 1.0, 0.5, 1.0, 0.5, 0.5, 0.0005, 0.0003];
  const DEFAULT_BOUNDS_SUP = [90.0, 90.0, 500.0, 10.0, 500.0, 10.0, 10.0, 0.0100, 0.0050];

  const el = {
    csv: document.getElementById('csv'),
    btnCalcular: document.getElementById('btnCalcular'),
    btnConfig: document.getElementById('btnConfig'),
    btnConfigFechar: document.getElementById('btnConfigFechar'),
    btnConfigRestore: document.getElementById('btnConfigRestore'),
    status: document.getElementById('status'),
    resultados: document.getElementById('resultados'),
    tabelaParams: document.getElementById('tabelaParams'),
    erroMae: document.getElementById('erroMae'),
    statsAviso: document.getElementById('statsAviso'),
    grafico: document.getElementById('grafico'),
    T_ini: document.getElementById('T_ini'),
    diametro: document.getElementById('diametro'),
    usarChute: document.getElementById('usarChute'),
    chuteContainer: document.getElementById('chuteContainer'),
    modalConfig: document.getElementById('modalConfig'),
    boundsBody: document.getElementById('boundsBody'),
    confianca: document.getElementById('confianca'),
    eps_rel: document.getElementById('eps_rel'),
    tabAnalise: document.getElementById('tabAnalise'),
    tabMonitoramento: document.getElementById('tabMonitoramento'),
    painelAnalise: document.getElementById('painelAnalise'),
    painelMonitoramento: document.getElementById('painelMonitoramento'),
    paramsAnaliseContainer: document.getElementById('paramsAnaliseContainer'),
    btnGerarCurva: document.getElementById('btnGerarCurva'),
    statusAnalise: document.getElementById('statusAnalise'),
    resultadosAnalise: document.getElementById('resultadosAnalise'),
    graficoAnalise: document.getElementById('graficoAnalise'),
    btnSalvarCSV: document.getElementById('btnSalvarCSV'),
    analise_T_ini: document.getElementById('analise_T_ini'),
    analise_diametro: document.getElementById('analise_diametro'),
    analise_C_cim: document.getElementById('analise_C_cim'),
    C_cim: document.getElementById('C_cim'),
    t_min: document.getElementById('t_min'),
    t_max: document.getElementById('t_max'),
    n_pontos: document.getElementById('n_pontos'),
  };

  let chartInstance = null;
  let chartAnaliseInstance = null;
  let lastCurvaData = null;

  function parseNum(val) {
    if (typeof val !== 'string') val = String(val || '');
    return parseFloat(val.replace(',', '.'));
  }

  function buildChuteInputs() {
    el.chuteContainer.innerHTML = '';
    NOMES_PARAMS.forEach((nome, i) => {
      const div = document.createElement('div');
      div.className = 'chute-item';
      div.innerHTML = `<label>${nome}</label><input type="number" step="any" data-chute-idx="${i}" value="${DEFAULT_CHUTE[i]}">`;
      el.chuteContainer.appendChild(div);
    });
  }

  function buildParamsAnaliseInputs() {
    if (!el.paramsAnaliseContainer) return;
    el.paramsAnaliseContainer.innerHTML = '';
    NOMES_PARAMS.forEach((nome, i) => {
      const div = document.createElement('div');
      div.className = 'chute-item';
      div.innerHTML = `<label>${nome}</label><input type="number" step="any" data-param-idx="${i}" value="${DEFAULT_CHUTE[i]}">`;
      el.paramsAnaliseContainer.appendChild(div);
    });
  }

  function getParamsAnalise() {
    const params = [];
    for (let i = 0; i < 9; i++) {
      const input = el.paramsAnaliseContainer.querySelector(`[data-param-idx="${i}"]`);
      params.push(input ? parseNum(input.value) : NaN);
    }
    return params.length === 9 && params.every(n => !Number.isNaN(n)) ? params : null;
  }

  function buildBoundsRows() {
    el.boundsBody.innerHTML = '';
    NOMES_PARAMS.forEach((nome, i) => {
      const tr = document.createElement('tr');
      tr.innerHTML = `
        <td>${nome}</td>
        <td><input type="number" step="any" data-bound-inf="${i}" value="${DEFAULT_BOUNDS_INF[i]}" autocomplete="off"></td>
        <td><input type="number" step="any" data-bound-sup="${i}" value="${DEFAULT_BOUNDS_SUP[i]}" autocomplete="off"></td>
      `;
      el.boundsBody.appendChild(tr);
    });
  }

  function getConfig() {
    const T_ini = parseNum(el.T_ini.value);
    const diametro = parseNum(el.diametro.value);
    const C_cim = parseNum(el.C_cim && el.C_cim.value);
    const confiancaPct = parseNum(el.confianca.value);
    const confianca = confiancaPct / 100;
    let eps_rel = parseNum(el.eps_rel.value);
    if (Number.isNaN(eps_rel)) eps_rel = 1e-4;

    const bounds_inf = [];
    const bounds_sup = [];
    el.boundsBody.querySelectorAll('tr').forEach((tr, i) => {
      const inf = tr.querySelector(`[data-bound-inf="${i}"]`);
      const sup = tr.querySelector(`[data-bound-sup="${i}"]`);
      if (inf) bounds_inf.push(parseNum(inf.value));
      if (sup) bounds_sup.push(parseNum(sup.value));
    });

    return {
      T_ini: Number.isNaN(T_ini) ? 25 : T_ini,
      diametro: Number.isNaN(diametro) ? 0.9 : diametro,
      C_cim: Number.isNaN(C_cim) ? 300 : C_cim,
      confianca: Number.isNaN(confianca) ? 0.95 : confianca,
      eps_rel,
      bounds_inf: bounds_inf.length === 9 ? bounds_inf : DEFAULT_BOUNDS_INF,
      bounds_sup: bounds_sup.length === 9 ? bounds_sup : DEFAULT_BOUNDS_SUP,
    };
  }

  function getChute() {
    if (!el.usarChute.checked) return undefined;
    const chute = [];
    for (let i = 0; i < 9; i++) {
      const input = el.chuteContainer.querySelector(`[data-chute-idx="${i}"]`);
      chute.push(input ? parseNum(input.value) : NaN);
    }
    return chute.length === 9 && chute.every(n => !Number.isNaN(n)) ? chute : undefined;
  }

  function setStatus(msg, isError) {
    el.status.textContent = msg || '';
    el.status.className = isError ? 'error' : (msg ? 'ok' : '');
  }

  function parseCSV(text) {
    const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
    if (lines.length < 2) return null;
    const rows = [];

    // Detecta o separador mais provável na primeira linha de dados válida
    const headerLine = lines[0];
    const dataLine = lines[1];
    let sep = ',';
    if (dataLine.includes(';')) sep = ';';
    else if (dataLine.includes('\t')) sep = '\t';

    for (let i = 0; i < lines.length; i++) {
      const parts = lines[i].split(sep);
      if (parts.length >= 2) {
        // Remove aspas eventuais e troca vírgula por ponto para o parseFloat
        const tStr = parts[0].replace(/['"]/g, '').replace(',', '.').trim();
        const TStr = parts[1].replace(/['"]/g, '').replace(',', '.').trim();
        const t = parseFloat(tStr);
        const T = parseFloat(TStr);
        // Só adiciona se for número válido (ignora cabeçalhos automaticamente)
        if (!isNaN(t) && !isNaN(T)) {
          rows.push({ t, T });
        }
      }
    }
    return rows.length ? rows : null;
  }

  function drawChart(res) {
    if (chartInstance) chartInstance.destroy();

    const dadosT = res.t_dados;
    const dadosTemp = res.T_dados;
    const tPlot = res.t_plot;
    const TPlot = res.T_plot;
    const ciLwr = res.CI_lwr;
    const ciUpr = res.CI_upr;

    const ctx = el.grafico.getContext('2d');
    chartInstance = new Chart(ctx, {
      type: 'line',
      data: {
        labels: tPlot.map(v => Number(v.toFixed(4))),
        datasets: [
          // 1. Limite Superior da Banda (invisível, serve de teto pro preenchimento)
          {
            label: 'IC 95% Superior',
            data: tPlot.map((t, i) => ({ x: t, y: ciUpr[i] })),
            borderColor: 'transparent',
            backgroundColor: 'transparent',
            pointRadius: 0,
            fill: false,
            tension: 0.1
          },
          // 2. Limite Inferior preenchendo até o dataset anterior (index 0)
          {
            label: 'Banda de Confiança ' + (res.confianca != null ? (res.confianca * 100).toFixed(0) : '95') + '%',
            data: tPlot.map((t, i) => ({ x: t, y: ciLwr[i] })),
            borderColor: 'transparent',
            backgroundColor: 'rgba(88, 166, 255, 0.2)',
            pointRadius: 0,
            fill: '-1',
            tension: 0.1
          },
          // 3. Curva Analítica Otimizada
          {
            label: 'Curva Analítica (Regressão)',
            data: tPlot.map((t, i) => ({ x: t, y: TPlot[i] })),
            borderColor: 'rgba(88, 166, 255, 1)',
            borderWidth: 2.5,
            pointRadius: 0,
            tension: 0.1,
          },
          // 4. Dados MEF Experimentais (Pontos)
          {
            label: 'Dados (Monitoramento)',
            data: dadosT.map((_, i) => ({ x: dadosT[i], y: dadosTemp[i] })),
            borderColor: 'rgb(248, 81, 73)',
            backgroundColor: 'transparent',
            borderWidth: 1,
            pointRadius: 2,
            pointBackgroundColor: 'rgb(248, 81, 73)',
            showLine: false, // Só exibe os pontos
          }
        ],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        interaction: { intersect: false, mode: 'index' },
        plugins: {
          legend: {
            position: 'top',
            labels: { filter: item => item.text !== 'IC 95% Superior' } // Esconde o teto da legenda
          },
        },
        scales: {
          x: { type: 'linear', title: { display: true, text: 'Tempo (h)' } },
          y: { title: { display: true, text: 'Temperatura (°C)' } },
        },
      },
    });
  }

  function drawChartCurva(t_plot, T_plot) {
    if (chartAnaliseInstance) chartAnaliseInstance.destroy();
    const ctx = el.graficoAnalise.getContext('2d');
    chartAnaliseInstance = new Chart(ctx, {
      type: 'line',
      data: {
        datasets: [{
          label: 'Temperatura (°C)',
          data: t_plot.map((t, i) => ({ x: t, y: T_plot[i] })),
          borderColor: 'rgba(88, 166, 255, 1)',
          backgroundColor: 'rgba(88, 166, 255, 0.1)',
          borderWidth: 2,
          pointRadius: 0,
          tension: 0.1,
        }],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: { legend: { position: 'top' } },
        scales: {
          x: { type: 'linear', title: { display: true, text: 'Tempo (h)' } },
          y: { title: { display: true, text: 'Temperatura (°C)' } },
        },
      },
    });
  }

  function showResultadosAnalise(res) {
    lastCurvaData = { t_plot: res.t_plot, T_plot: res.T_plot };
    drawChartCurva(res.t_plot, res.T_plot);
    el.resultadosAnalise.classList.add('visible');
  }

  function runCurva() {
    const params = getParamsAnalise();
    if (!params) {
      el.statusAnalise.textContent = 'Preencha os 9 parâmetros.';
      el.statusAnalise.className = 'error';
      return;
    }
    const T_ini = parseFloat(el.analise_T_ini && el.analise_T_ini.value) || 25;
    const diametro = parseFloat(el.analise_diametro && el.analise_diametro.value) || 0.9;
    const C_cim = parseFloat(el.analise_C_cim && el.analise_C_cim.value) || 300;
    const tMin = parseFloat(el.t_min && el.t_min.value) || 0.1;
    const tMax = parseFloat(el.t_max && el.t_max.value) || 100;
    const n = Math.max(50, Math.min(2000, parseInt(el.n_pontos && el.n_pontos.value, 10) || 300));
    const tempos = [];
    for (let i = 0; i < n; i++) tempos.push(tMin + (tMax - tMin) * i / (n - 1));
    el.statusAnalise.textContent = 'Gerando curva...';
    el.statusAnalise.className = 'ok';
    el.btnGerarCurva.disabled = true;
    fetch(API_CURVA, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        params,
        config: { T_ini, diametro, C_cim },
        tempos,
      }),
    })
      .then(async r => {
        const text = await r.text();
        try {
          const data = JSON.parse(text);
          if (!r.ok) throw new Error(data.error || 'Erro na API');
          return data;
        } catch (e) {
          console.error("Resposta do servidor:", text);
          if (text.includes('504')) throw new Error('Tempo esgotado (Timeout). O servidor do Render é lento demais para este cálculo.');
          throw new Error('Erro do Servidor (Não JSON). Verifique os logs do Render.');
        }
      })
      .then(showResultadosAnalise)
      .then(() => {
        el.statusAnalise.textContent = 'Pronto.';
        el.btnGerarCurva.disabled = false;
      })
      .catch(err => {
        el.statusAnalise.textContent = err.message || 'Erro.';
        el.statusAnalise.className = 'error';
        el.btnGerarCurva.disabled = false;
      });
  }

  function saveCurvaCSV() {
    if (!lastCurvaData) return;
    const { t_plot, T_plot } = lastCurvaData;
    const lines = ['tempo_h,temperatura_C'];
    t_plot.forEach((t, i) => lines.push(t + ',' + T_plot[i]));
    const blob = new Blob([lines.join('\n')], { type: 'text/csv;charset=utf-8' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'curva_tubulao.csv';
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function showResultados(res) {
    el.erroMae.textContent = 'RMSE / MAE do ajuste: ' + (res.erro_mae != null ? res.erro_mae.toFixed(4) : '—') + ' °C';
    if (res.stats_aviso && el.statsAviso) {
      el.statsAviso.textContent = res.stats_aviso + (res.cond_FtF != null ? ' [cond(F\'F) = ' + Number(res.cond_FtF).toExponential(2) + ']' : '');
      el.statsAviso.style.display = 'block';
      el.statsAviso.style.color = 'var(--error)';
    } else if (res.cond_FtF != null && el.statsAviso) {
      el.statsAviso.textContent = 'cond(F\'F) = ' + Number(res.cond_FtF).toExponential(2);
      el.statsAviso.style.display = 'block';
      el.statsAviso.style.color = 'var(--muted)';
    } else if (el.statsAviso) {
      el.statsAviso.style.display = 'none';
    }
    const pct = res.confianca != null ? (res.confianca * 100).toFixed(0) : '95';
    const thIcInf = document.getElementById('thIcInf');
    const thIcSup = document.getElementById('thIcSup');
    if (thIcInf) thIcInf.textContent = 'IC Inf (' + pct + '%)';
    if (thIcSup) thIcSup.textContent = 'IC Sup (' + pct + '%)';

    // Preenche a tabela
    const params = res.parametros || [];
    el.tabelaParams.innerHTML = params.map(p => `
      <tr>
        <td>${p.nome}</td>
        <td>${p.estimado.toFixed(4)}</td>
        <td>${p.se.toFixed(4)}</td>
        <td>${p.ic_inf.toFixed(4)}</td>
        <td>${p.ic_sup.toFixed(4)}</td>
        <td style="color: ${p.cv > 30 ? 'var(--error)' : 'inherit'}">${p.cv.toFixed(1)}%</td>
      </tr>
    `).join('');

    drawChart(res);
    el.resultados.classList.add('visible');
  }

  function runOtimizacao(tempos, temperaturas) {
    setStatus('Calculando (Jacobiano 9D)... Por favor aguarde, pode levar até 1 minuto.', false);
    el.btnCalcular.disabled = true;

    const config = getConfig();
    const chute = getChute();
    const payload = { tempos, temperaturas, config };
    if (chute) payload.chute = chute;

    fetch(API_OTIMIZAR, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    })
      .then(async r => {
        const text = await r.text();
        try {
          const data = JSON.parse(text);
          if (!r.ok) throw new Error(data.error || 'Erro na API');
          data.t_dados = tempos;
          data.T_dados = temperaturas;
          return data;
        } catch (e) {
          console.error("Resposta do servidor:", text);
          if (text.includes('504')) throw new Error('Tempo esgotado (Timeout). Render é muito lento para essa regressão 9D.');
          throw new Error('Erro do Servidor (Não JSON). Verifique os logs no painel do Render.');
        }
      })
      .then(showResultados)
      .then(() => {
        setStatus('Pronto.');
        el.btnCalcular.disabled = false;
      })
      .catch(err => {
        setStatus(err.message || 'Erro ao chamar a API.', true);
        el.btnCalcular.disabled = false;
      });
  }

  el.usarChute.addEventListener('change', () => {
    el.chuteContainer.style.display = el.usarChute.checked ? 'grid' : 'none';
  });

  el.btnConfig.addEventListener('click', () => {
    el.modalConfig.classList.add('visible');
  });
  el.btnConfigFechar.addEventListener('click', () => {
    el.modalConfig.classList.remove('visible');
  });
  el.modalConfig.addEventListener('click', e => {
    if (e.target === el.modalConfig) el.modalConfig.classList.remove('visible');
  });

  el.btnConfigRestore.addEventListener('click', () => {
    el.confianca.value = 95;
    el.eps_rel.value = '1e-4';
    el.boundsBody.querySelectorAll('tr').forEach((tr, i) => {
      const inf = tr.querySelector(`[data-bound-inf="${i}"]`);
      const sup = tr.querySelector(`[data-bound-sup="${i}"]`);
      if (inf) inf.value = DEFAULT_BOUNDS_INF[i];
      if (sup) sup.value = DEFAULT_BOUNDS_SUP[i];
    });
  });

  buildChuteInputs();
  buildBoundsRows();
  buildParamsAnaliseInputs();

  el.tabAnalise.addEventListener('click', () => {
    el.tabAnalise.classList.add('active');
    el.tabMonitoramento.classList.remove('active');
    el.painelAnalise.classList.add('visible');
    el.painelMonitoramento.classList.remove('visible');
  });
  el.tabMonitoramento.addEventListener('click', () => {
    el.tabMonitoramento.classList.add('active');
    el.tabAnalise.classList.remove('active');
    el.painelMonitoramento.classList.add('visible');
    el.painelAnalise.classList.remove('visible');
  });

  el.btnGerarCurva.addEventListener('click', runCurva);
  el.btnSalvarCSV.addEventListener('click', saveCurvaCSV);

  el.btnCalcular.addEventListener('click', () => {
    const file = el.csv.files && el.csv.files[0];
    if (!file) {
      setStatus('Selecione um arquivo CSV.', true);
      return;
    }
    const reader = new FileReader();
    reader.onload = e => {
      const rows = parseCSV(e.target.result);
      if (!rows) {
        setStatus('CSV inválido ou sem dados numéricos.', true);
        return;
      }
      const tempos = rows.map(r => r.t);
      const temperaturas = rows.map(r => r.T);
      runOtimizacao(tempos, temperaturas);
    };
    reader.readAsText(file, 'UTF-8');
  });
})();