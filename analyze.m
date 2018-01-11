pkg load optim % leasqr
pkg load tmvs % filters, mapl

models = {'none', 'hw', 'cs'};
ps = {'0.0125', '0.025', '0.0375', '0.05', '0.0625', ...
  '0.075', '0.0875', '0.1', '0.1125', '0.125', ...
  '0.15', '0.175', '0.2', '0.225', '0.25', ...
  '0.3', '0.35', '0.4', '0.45', '0.5', ...
  '0.6', '0.7', '0.8', '0.9', '1.0'};

if (!exist ('runs', 'var'))
  runs = struct ('model', {}, 'p', {}, ...
    't', {}, 'e', {}, 'q', {}, 'w', {}, ...
    'mu', {}, 'v', {}, 'f', {});

  for imodel = 1 : length (models)
    model = models{imodel};

    for ip = 1 : length (ps)
      p = ps{ip};

      [t, e, q, w, mudir, mufb, vx, vy, fx, fy] = textread ( ...
        sprintf ('run/%s-%s/est.data', model, p));

      runs(length (runs) + 1) = struct ('model', model, 'p', str2num (p), ...
        't', t - t(1), 'e', e, 'q', q, 'w', w, ...
        'mu', [mudir, mufb], 'v', [vx, vy], 'f', [fx, fy]);
    end
  end
end

fexp = @(x0, y0) @(x, q) ifelse (x < x0, (y0 / x0) * x, ...
  q(1) + (y0 - q(1)) * exp (-q(2) * (x - x0)));

if (!exist ('fits', 'var'))
  fits = struct ('run', {}, 't', {}, 'mu', {}, 'q', {}, 'dq', {});

  for imodel = 1 : length (models)
    model = models{imodel};

    fruns = filters (@(run) strcmp (run.model, model), runs);

    for irun = 1 : length (fruns)
      run = fruns(irun);

      x = run.t;
      y = run.mu(:, 1);
      y0 = max (y);
      i0 = find (y == y0)(1);
      x0 = x(i0);
      xcut = x(i0 : end);
      ycut = y(i0 : end);
      [~, q] = leasqr (xcut, ycut, [0.25, 1.5e+3], fexp (x0, y0));

      fits(length (fits) + 1) = struct ('run', run, ...
        't', xcut, 'mu', ycut, 'q', q, 'dq', sqrt ( ...
        sum ((ycut - fexp (x0, y0) (xcut, q)) .^ 2) / (length (ycut) - 1)));
    end
  end
end

flin = @(x0, y0) @(x, q) (q(1) + y0) + q(2) * (x - x0);

if (!exist ('exts', 'var'))
  exts = struct ('run', {}, 't', {}, 'w', {}, 'q', {}, 'dq', {});

  for imodel = 1 : length (models)
    model = models{imodel};

    fruns = filters (@(run) strcmp (run.model, model), runs);

    for irun = 1 : length (fruns)
      run = fruns(irun);

      x = run.t;
      y = run.w;
      i0 = x >= mean ([(min (x)), (max (x))]);
      xcut = x(i0);
      ycut = y(i0);
      x0 = xcut(1);
      y0 = ycut(1);
      [~, q] = leasqr (xcut, ycut, [0.0, 1.0e+3], flin (x0, y0));

      exts(length (exts) + 1) = struct ('run', run, ...
        't', xcut, 'w', ycut, 'q', q, 'dq', sqrt ( ...
        sum ((ycut - flin (x0, y0) (xcut, q)) .^ 2) / (length (ycut) - 1)));
    end
  end
end

fps = {'0.15', '0.6'};

for imodel = 1 : length (models)
  model = models{imodel};

  ffits = filters (@(fit) strcmp (fit.run.model, model), fits);

  figure (imodel);
  clf ();
  hold ('on');
  h = [];

  for ifit = 1 : length (ffits)
    fit = ffits(ifit);

    if (find (fit.run.p == cell2mat (mapl (@str2num, fps))))
      c = [(interp1 ([0.0, 0.5, 1.0], [0.0, 0.5, 1.0], fit.run.p)), ...
        (interp1 ([0.0, 0.5, 1.0], [0.0, 0.0, 0.0], fit.run.p)), ...
        (interp1 ([0.0, 0.5, 1.0], [1.0, 0.5, 0.0], fit.run.p))];

      plot (fit.run.t, fit.run.mu, ...
        'color', c,  'linewidth', 1);

      h(length (h) + 1) = plot (fit.run.t, ...
        fexp (fit.t(1), fit.mu(1)) (fit.run.t, fit.q), ...
        'color', c, 'linewidth', 2);

      for isgn = 1 : 2
        sgn = [-1.0, 1.0](isgn);

        plot (fit.t, ...
          fexp (fit.t(1), fit.mu(1)) (fit.t, fit.q) + sgn * fit.dq, ...
          'color', c, 'linewidth', 1);
      end
    end
  end

  title (model);
  legend (h, mapl (@(p) ['p = ', p], fps), 'location', 'northeast');
  xlabel ('t');
  ylabel ('\mu');
  axis ([0.0e-3, 50.0e-3, 0.0, 2.0]);
  hold ('off');
end

figure (length (models) + 1);
clf ();
hold ('on');
h = [];

for imodel = 1 : length (models)
  model = models{imodel};

  ffits = filters (@(fit) strcmp (fit.run.model, model), fits);

  c = 0.5 * [imodel == 1, imodel == 2, imodel == 3];

  h(length (h) + 1) = plot ([[ffits.run].p], [ffits.q](1, :), ...
      'color', c, 'linewidth', 2);

  for isgn = 1 : 2
    sgn = [-1.0, 1.0](isgn);

    plot ([[ffits.run].p], [ffits.q](1, :) + sgn * [ffits.dq], ...
      'color', c, 'linewidth', 1);
  end
end

legend (h, models, 'location', 'northeast');
xlabel ('p / p^{crit}');
ylabel ('\mu');
axis ([0.0, 0.6, 0.0, 1.0]);
hold ('off');

fps = {'0.0625', '0.125', '0.25', '0.5'};

for ip = 1 : length (fps)
  p = fps{ip};

  fexts = filters (@(ext) ext.run.p == str2num (p), exts);

  figure (length (models) + 1 + ip);
  clf ();
  hold ('on');
  h = [];

  for iext = 1 : length (fexts)
    ext = fexts(iext);

    model = ext.run.model;
    imodel = find (strcmp (models, model));

    c = 0.5 * [imodel == 1, imodel == 2, imodel == 3];

    plot (ext.run.t, ext.run.w, ...
      'color', c, 'linewidth', 1);

    h(length (h) + 1) = plot (ext.t, ...
      flin (ext.t(1), ext.w(1)) (ext.t, ext.q), ...
      'color', c, 'linewidth', 2);

    for isgn = 1 : 2
      sgn = [-1.0, 1.0](isgn);

      plot (ext.t, ...
        flin (ext.t(1), ext.w(1)) (ext.t, ext.q) + sgn * ext.dq, ...
        'color', c, 'linewidth', 1);
    end
  end

  title (['p = ', p]);
  legend (h, models, 'location', 'northwest');
  xlabel ('t');
  ylabel ('W');
  axis ([0.0e-3, 50.0e-3]);
  hold ('off');
end
