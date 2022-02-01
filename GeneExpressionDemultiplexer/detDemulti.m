% ODEs describing the demultiplexer
% y(1) = Activator1
% y(2) = Repressor
% y(3) = target mRNA1
% y(4) = target Protein1
% y(5) = Activator2
% y(6) = target mRNA2
% y(7) = target Protein2


function dy = detDemulti(t,y,p,p2,TFtot,TFtot2,Reptot,I,kdegProt)
dy = zeros(7,1);
dy(1) = (TFtot - y(1)) * p(1) * I - y(1) * p(2);
dy(2) = (Reptot - y(2)) * p(3) * I - y(2) * p(4);
dy(3) = p(5) + p(6) * y(1) ^ p(8) / (y(1) ^ p(8) + p(7) ^ p(8))  * (1 /(1 + (y(2)/p(9))^p(10))) - y(3) * p(11);
dy(4) = p(12) * y(3) - kdegProt * y(4);
dy(5) = (TFtot2 - y(5)) * p2(1) * I - y(5) * p2(2);
dy(6) = p2(3) + p2(4) * y(5) ^ p2(6) / (y(5) ^ p2(6) + p2(5) ^ p2(6)) - y(6) * p2(7);
dy(7) = p2(8) * y(6) - kdegProt * y(7);


%p are the parameters of the diamond-IFFL

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

%p2 are the parameters of the second gene expression system

%p2(1) -> on rate Act
%p2(2) -> off rate Act
%p2(3) -> basal transcription
%p2(4) -> max transcription
%p2(5) -> MM constant Act
%p2(6) -> Hill coeff Act
%p(7) -> mRNA degradation rate
%p(8) -> translation rate / mRNA

end

