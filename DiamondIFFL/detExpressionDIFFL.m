% ODEs describing the diamond-IFFL
% y(1) = Activator
% y(2) = Repressor
% y(3) = target mRNA
% y(4) = target Protein

function dy = detExpressionDIFFL(t,y,p,TFtot,Reptot,I,kdegProt)
dy = zeros(4,1);
dy(1) = (TFtot - y(1)) * p(1) * I - y(1) * p(2);
dy(2) = (Reptot - y(2)) * p(3) * I - y(2) * p(4);
dy(3) = p(5) + p(6) * y(1) ^ p(8) / (y(1) ^ p(8) + p(7) ^ p(8))  * (1 /(1 + (y(2)/p(9))^p(10))) - y(3) * p(11);
dy(4) = p(12) * y(3) - kdegProt * y(4);

%p(1) -> on rate Act
%p(2) -> off rate Act
%p(3) -> on rate Rep
%p(4) -> off rate Rep
%p(5) -> basal transcription
%p(6) -> max transcription
%p(7) -> MM constant Act
%p(8) -> Hill coeff Act
%p(9) -> MM constant Rep
%p(10) -> Hill coeff Rep
%p(11) -> mRNA degradation rate
%p(12) -> translation rate / mRNA
%kdegProt -> Fluorescent protein degradation rate

end

