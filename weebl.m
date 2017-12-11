[x, y] = textread ('trash.data');
[~, p] = leasqr (x, y, [1.6e-3, 6.4], @(x, p) wblpdf (x, num2cell (p){:}));
plot (x, y, x, wblpdf (x, num2cell (p){:}));
