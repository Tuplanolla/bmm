pkg load data-smoothing
pkg load optim
pkg load tmvs

if (!exist ('runs', 'var'))
  models = {'none', 'hw', 'cs'};
  ps = {'0.01', ...
    '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'};

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

someruns = filters (@(run) strcmp (run.model, 'hw'), runs);

plim = [];
mulim = [];

figure (1);
clf ();
xlabel ('t');
ylabel ('\mu');
axis ([0.0e-3, 50.0e-3, 0.0, 2.0]);
hold ('on');

for irun = 1 : length (someruns)
  run = someruns(irun);

  if (run.p >= 0.01)
    x = run.t;
    y = run.mu(:, 2);
    y0 = max (y);
    i0 = find (y == y0)(1);
    x0 = x(i0);
    xcut = x(i0 : end);
    ycut = y(i0 : end);
    f = @(x, q) q(1) + (y0 - q(1)) * exp (-q(2) * (x - x0));
    [~, q] = leasqr (xcut, ycut, [0.25, 1.5e+3], f);

    plim(irun) = run.p;
    mulim(irun) = q(1);

    c = [run.p, 0.0, 1.0 - run.p];
    plot (xcut, f (xcut, q), 'color', c);
    % plot (x, y, 'color', c);
    plot (x, regdatasmooth (x, y, 'lambda', 1.0e-6), 'color', c);
  end
end

hold ('off');

figure (2);
clf ();
xlabel ('p');
ylabel ('\mu');
axis ([0.0, 1.0, 0.0, 2.0]);
hold ('on');

plot (plim, mulim);

hold ('off');
