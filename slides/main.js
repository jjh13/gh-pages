$(function(){
	ko.bindingHandlers.slider = {
	init: function (element, valueAccessor, allBindingsAccessor) {
		var options = allBindingsAccessor().sliderOptions || {};
		$(element).slider(options);
		ko.utils.registerEventHandler(element, "slidechange", function (event, ui) {
			var observable = valueAccessor();
			observable(ui.value);
		});
		ko.utils.domNodeDisposal.addDisposeCallback(element, function () {
			$(element).slider("destroy");
		});
		ko.utils.registerEventHandler(element, "slide", function (event, ui) {
			var observable = valueAccessor();
			observable(ui.value);
		});
	},
	update: function (element, valueAccessor) {
			var value = ko.utils.unwrapObservable(valueAccessor());
			if (isNaN(value)) value = 0;
			$(element).slider("value", value);
		}
	};

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
		time: 100,
	};

	function heaviside(t) {
		if (t > 0) return 1.;
		return 0;
	};

	function intrahost(t,x) {
		var T = x[0], I1 = x[1], I2 = x[2], I3 = x[3],	I4 = x[4], I5 = x[5], I6 = x[6], M1 = x[7],
			M2 = x[8], L1 = x[9], L2 = x[10], D = x[11], V = x[12];

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
	};

	/**
	 * Graphing
	 */

	var palette = new Rickshaw.Color.Palette();
	var graph_t = new Rickshaw.Graph( {
    	    element: document.querySelector("#chart_t"),
        	width: 580,
        	height: 250,
        	renderer: 'line',
        	series: [
        		{name: "T-cells" , color: palette.color(), data: [] }, 
        		{name: "DNA+ cells", color: palette.color(), data: [] }, ]
	});

	// Bang off one more to make things prettier
	palette.color();


	var graph_v = new Rickshaw.Graph( {
    	    element: document.querySelector("#chart_v"),
        	width: 580,
        	height: 250,
        	renderer: 'line',
        	series: [
        		{name:"Low", color: palette.color(), data: [] }, 
        		{name:"Medium",color: palette.color(), data: [] }, 
        		{name:"High",color: palette.color(), data: [] }, 
        		{name:"Extra",color: palette.color(), data: [] },]
	});

	var x_unit = {
		name: 'day',
			seconds: 86400 * 60,
			formatter: function(d) { return Math.floor(d/8.64e7); }
	}
	new Rickshaw.Graph.Axis.Time( { graph: graph_t, timeUnit: x_unit });
	new Rickshaw.Graph.Axis.Time( { graph: graph_v, timeUnit: x_unit });

	new Rickshaw.Graph.Axis.Y( {
		graph: graph_t,
		orientation: 'left',
		tickFormat: Rickshaw.Fixtures.Number.formatKMBT,
		element: document.getElementById('y_axis_t'),
	} );

	new Rickshaw.Graph.Axis.Y( {
		graph: graph_v,
		orientation: 'left',
		tickFormat: Rickshaw.Fixtures.Number.formatKMBT,
		element: document.getElementById('y_axis_v'),
	});

	new Rickshaw.Graph.Legend({ element: document.querySelector('#legend_t'), graph: graph_t});
	new Rickshaw.Graph.Legend({ element: document.querySelector('#legend_v'), graph: graph_v});
	new Rickshaw.Graph.HoverDetail({graph: graph_t});
	new Rickshaw.Graph.HoverDetail({graph: graph_v});

	function prep_output(x, y) {
		var data = new Array(x.length);

		for(var i = 0; i < x.length; i++)
			data[i] = {x:x[i]*86400, y:y[i]};
	
		return data;
	};

	var ViewModel = function() {
		var self = this;
		var add = numeric.add;
		
		self.effectiveness = ko.observable(0.85);
		self.tttreatment = ko.observable(365);
		
		self.net = ko.computed(function() {
			p.time = self.tttreatment();
			p.epsilon = self.effectiveness();

			var sol = numeric.dopri(0, 2000 ,[p.lambda/p.delta_t, 0,0,0,0,0,0,0,0,0,0,0, 1],intrahost,1e-3,10000);
			var y = numeric.transpose(sol.y);

			graph_t.series[0].data = prep_output(sol.x, y[0]);

			graph_t.series[1].data = prep_output(sol.x, 
				add(add(add(y[2], y[3]), add(y[4], y[5])), 
				add(add(add(y[6], y[7]), add(y[8], y[9])), add(y[10], y[11]))));

			graph_v.series[0].data = prep_output(sol.x, add(add(y[1], y[3]), y[9]));
			graph_v.series[1].data = prep_output(sol.x, add(add(y[4], y[10]), y[7]));
			graph_v.series[2].data = prep_output(sol.x, add(add(y[5], y[6]), y[8]));
			graph_v.series[3].data = prep_output(sol.x, add(y[6], y[8]));
			graph_t.render();
			graph_v.render();
			return null;
		});
	}

	ko.applyBindings(new ViewModel());
});