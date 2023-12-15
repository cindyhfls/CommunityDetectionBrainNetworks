function input_check_HSB(clu)
    assert(sum(all(clu==0))==0,'some columns are full of 0s, considering replacing it with 1')
    assert(sum(isnan(clu(:)))==0,'input contains NaN');  
end