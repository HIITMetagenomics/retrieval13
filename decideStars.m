function stars = decideStars(pvalue)

stars = 0;
if pvalue < 0.05 && pvalue >= 0.01
    stars = 1;
end
if pvalue < 0.01 && pvalue >= 0.001
    stars = 2;
end
if pvalue < 0.001 && pvalue >= 0.0
    stars = 3;
end