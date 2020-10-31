function ind = box_co(a,b,dim)
dim2 = dim / 2;
if a <= dim2
    if b <= dim2
        ind = 3;
    else
        ind = 4;
        b = b - dim2;
    end
else
    a = a - dim2;
    if b <= dim2
        ind = 1;
    else
        ind = 2;
        b = b - dim2;
    end
end
if dim ~= 2
    ind = [ind, box_co(a, b, dim2)];
end
 

