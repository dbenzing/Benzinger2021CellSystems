%Function to calculate response of diamond-IFFL to PWM (defined by period, width, light intensity (I)) 
%for given parameter set (p, TFtot, Reptot, kdegProt) and experimental timespan (tspan) 

function PWMss = PWMsteady(p,initial,TFtot,Reptot,I, period, width,kdegProt,tspan)

	pulsenumber = tspan(end) / period;

	for numb = 1 : pulsenumber
		tspan = [0 width];
		[T,Y] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt), tspan, initial);
		initial = Y(end,:);
			
		tspan = [0 period-width];
		[T,Y] = ode23s(@(t,y) detExpressionDIFFL(t,y,p,TFtot,Reptot,0,kdegProt), tspan, initial);
		initial = Y(end,:);
	end
	
	PWMss = Y(end,4);
end


