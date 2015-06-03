var p = {
    lambda: 68.4,
    epsilon: 0.85,
    beta: 1.4e-6,
    delta_t: 0.017,
    delta_d: 8.3e-5,

    delta_l1: 8.3e-5,
    delta_l2: 8.3e-5,

    delta_i1: 0.017,
    delta_i2: 0.017,
    delta_i3: 0.017,
    delta_i4: 0.017,
    delta_i5: 0.017,
    delta_i6: 1.0,

    delta_m1: 0.017,
    delta_m2: 1.0,

    c: 23.0,
    f_d: 0.141,
    f_l: 0.004,
    f_m: 0.25,
    gamma1: 3,
    gamma2: 6,
    gamma3: 12,
    gamma4: 12,
    gamma5: 12,
    sigma1: 0.079,
    sigma2: 1.0,
    k1: 0.103,
    k2: 1.0,
    alpha: 2.7e-3,
    n: 2.14e4,
    q: 2.14e4,
    time: 80,
};

heaviside = function(t){
	if (t > 0) return 1.;
	return 0;
}

f = function(t,x) {
	var T = x[0],
		I1 = x[1],
		I2 = x[2],
		I3 = x[3],
		I4 = x[4],
		I5 = x[5],
		I6 = x[6],
		M1 = x[7],
		M2 = x[8],
		L1 = x[9],
		L2 = x[10],
		D = x[11],
		V = x[12];

	return [
  			// dT/dt
  			p.lambda - (1 - p.epsilon*heaviside(t-p.time))*p.beta*V*T - p.delta_t*T,
  			
  			// dI1/dt
  			(1 - p.epsilon*heaviside(t-p.time))*p.beta*V*T - (p.gamma1 + p.delta_i1)*I1,
  			
  			// dI2/dt
  			p.gamma1*I1 - (p.gamma2 + p.delta_i2)*I2,
  			
  			// dI3/dt
  			(1-p.f_d)*p.gamma2*I2 - (p.gamma3 + p.delta_i3)*I3,
  			
  			// dI4/dt
  			(1-p.f_l)*p.gamma3*I3 + p.alpha*L2 - (p.gamma4 + p.delta_i4)*I4,

  			// dI5/dt
  			(1-p.f_m)*p.gamma4*I4 - (p.gamma5 + p.delta_i5)*I5,

  			// dI6/dt
  			p.gamma5*I5 - p.delta_i6*I6,

  			// dM1/dt
  			p.f_m*p.gamma4*I4 + p.k2*M2 - (p.k1 + p.delta_m1)*M1,

  			// dM2/dt
  			p.k1*M1 - (p.k2 + p.delta_m2)*M2,

  			// dL1/dt
  			p.f_l*p.gamma3*I3  + p.sigma2*L2 - (p.sigma1 + p.delta_l1)*L1,

  			// dL2/dt
  			p.sigma1*L1 - (p.alpha + p.sigma2 + p.delta_l2)*L2,

  			// dD/dt 
  			p.f_d*p.gamma2*I2 - p.delta_d*D,

  			// dV/dt
  			p.n*p.delta_i6*I6 + p.q*M2 - p.c*V ];
}

[p.lambda/p.delta_t, 0,0,0,0,0,0,0,0,0,0,0, 1]
