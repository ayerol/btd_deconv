function cmap = custom_colormap(mn,mx)

numColors = 100;
c = colormap(colorcet('CBTD1')); 
start_color = c(end,:); 
end_color = [.1 .6 .1];
pos_count = round(numColors*abs(mn)/(abs(mn)+mx))-2;
neg_count = numColors-pos_count;
r = linspace(start_color(1),1,pos_count);
r2 = linspace(1,end_color(1),neg_count+1); r2(1) = [];
r = [r r2]';
g = linspace(start_color(2),1,pos_count);
g2 = linspace(1,end_color(2),neg_count+1); g2(1) = [];
g = [g g2]';
b = linspace(start_color(3),1,pos_count);
b2 = linspace(1,end_color(3),neg_count+1); b2(1) = [];
b = [b b2]';
cmap = [r g b];

end