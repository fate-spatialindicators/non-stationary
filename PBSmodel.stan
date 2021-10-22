
	data 
		{
		int<lower=1> n;
		int<lower=1> num_knots;
		int<lower=1> degree;
		vector[n] y;
		matrix[n, num_knots+degree] B;
		matrix[num_knots+degree, num_knots+degree] Q;
		}
		
	parameters 
		{
		vector[num_knots+degree] theta;
		real<lower=0> sigma2;
		}
		
	transformed parameters
		{
		vector[n] mu;
		mu = B * theta;
		}
		
	model 
		{
		sigma2 ~ inv_chi_square(1);
		target += -0.5 * quad_form(Q, theta)/sigma2;
		y ~ normal(mu, sqrt(sigma2));
		}
		
	generated quantities
		{
		vector[n] yhat;
		yhat = B * theta;
		}
	