function norm_sqr = myNormSqr(k,whichBasis)
     switch whichBasis
     case 'Legendre'
        norm_sqr = 2./(2.*(k-1)+1);
     otherwise
     end
