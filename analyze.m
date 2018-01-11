pkg load optim
pkg load tmvs

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
        sprintf ('%s-%s/est.data', model, p));

      runs(length (runs) + 1) = struct ('model', model, 'p', str2num (p), ...
        't', t - t(1), 'e', e, 'q', q, 'w', w, ...
        'mu', [mudir, mufb], 'v', [vx, vy], 'f', [fx, fy]);
    end
  end
end

f = @(x0, y0) @(x, q) q(1) + (y0 - q(1)) * exp (-q(2) * (x - x0));

if (!exist ('fits', 'var'))
  fits = struct ('run', {}, 't', {}, 'mu', {}, 'q', {}, 'dq', {});

  for imodel = 1 : length (models)
    model = models{imodel};

    fruns = filters (@(run) strcmp (run.model, model), runs);

    figure (imodel);
    clf ();
    xlabel ('t');
    ylabel ('\mu');
    axis ([0.0e-3, 50.0e-3, 0.0, 2.0]);
    hold ('on');

    for irun = 1 : length (fruns)
      run = fruns(irun);

      x = run.t;
      y = run.mu(:, 1);
      y0 = max (y);
      i0 = find (y == y0)(1);
      x0 = x(i0);
      xcut = x(i0 : end);
      ycut = y(i0 : end);
      [~, q] = leasqr (xcut, ycut, [0.25, 1.5e+3], f (x0, y0));

      fits(length (fits) + 1) = struct ('run', run, ...
        't', xcut, 'mu', ycut, 'q', q, 'dq', sqrt ( ...
        sum ((ycut - f (x0, y0) (xcut, q)) .^ 2) / (length (ycut) - 1)));
    end
  end
end

for imodel = 1 : length (models)
  model = models{imodel};

  ffits = filters (@(fit) strcmp (fit.run.model, model), fits);

  figure (imodel);
  clf ();
  title (model);
  xlabel ('t');
  ylabel ('\mu');
  axis ([0.0e-3, 50.0e-3, 0.0, 2.0]);
  hold ('on');

  for ifit = 1 : length (ffits)
    fit = ffits(ifit);

    c = [fit.run.p, 0.0, 1.0 - fit.run.p];

    plot (fit.t, f (fit.t(1), fit.mu(1)) (fit.t, fit.q), ...
      'color', c, 'linewidth', 2);

    plot (fit.run.t, fit.run.mu, ...
      'color', c, 'linewidth', 1);
  end

  hold ('off');
end

figure (length (models) + 1);
clf ();
% legend (models);
xlabel ('p / p^{crit}');
ylabel ('\mu');
axis ([0.0, 1.0, 0.0, 1.0]);
hold ('on');

for imodel = 1 : length (models)
  model = models{imodel};

  ffits = filters (@(fit) strcmp (fit.run.model, model), fits);

  c = 0.5 * [imodel == 1, imodel == 2, imodel == 3];

  plot ([[ffits.run].p], [ffits.q](1, :), ...
      'color', c, 'linewidth', 2);

  for isgn = 1 : 2
    sgn = [-1.0, 1.0](isgn);

    plot ([[ffits.run].p], [ffits.q](1, :) + sgn * [ffits.dq], ...
      'color', c, 'linewidth', 1);
  end
end

hold ('off');
