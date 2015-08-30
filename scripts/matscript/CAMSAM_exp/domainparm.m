function [nsubdomain_x,nsubdomain_y,nx_gl,ny_gl] = domainparm(domain_x,domain_y,subdomain_x,subdomain_y,dx,dy)

if (rem(domain_x,dx)==0 && rem(domain_y,dy)==0)
   nx_gl = domain_x/dx;
   ny_gl = domain_y/dx;
   if (any(factor(nx_gl)>5) && any(factor(ny_gl)>5)) 
      error('nx_gl or ny_gl must be a factor of 2,3, or 5')
   end
   if (rem(domain_x,subdomain_x)==0 && rem(domain_y,subdomain_y)==0)
      nsubdomain_x = domain_x/subdomain_x;
      nsubdomain_y = domain_y/subdomain_y;
   else
      error('can not divide domain into integer number of subdomains, reset domain size')
   end
else
   error('cannot divide domain into integers by dx & dy, reset domain size')
end

