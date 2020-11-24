function [q] = normal_quantile(p)

    %https://en.wikipedia.org/wiki/Normal_distribution#Quantile_function

	pivot = 0.09;

	n_length = length(p);

	q = zeros(1,n_length);

	F = sqrt( -2*log(pivot) - log( -2*log(pivot) ) - log(2*pi) );
    
	for i = 1:n_length 
        
        if (p(i) < pivot)
			q(i) = - sqrt( -2*log(p(i)) - log( -2*log(p(i)) ) - log(2*pi) ) ;
        elseif (p(i) > 1 - pivot)
			q(i) = sqrt( -2*log(1-p(i)) - log( -2*log(1-p(i)) ) - log(2*pi) ) ;
        else
			a1 = 1/(2*pivot - 1);
			a2 = 1/(1-2*pivot);
			b1 = (1-pivot)/(1-2*pivot);
			b2 = 1/(2 - (1/pivot));
			
			q(i) = (a1*p(i)+b1)*(-F) + (a2*p(i)+b2)*F;
            
        end
    end
    
end