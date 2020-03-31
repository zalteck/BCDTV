function [zminus, eminus] = computingEsZs(YT,CT,M)

ns = size(CT,1);
nc = size(YT,1);
tamm = size(YT,2);
eminus = zeros(tamm, nc, ns);
zminus = zeros(tamm, ns);

%tic

YT_MCT = YT - M * CT;
for s=1:ns
    aux =  YT_MCT + M(:,s) * CT(s,:);
    eminus(:,:,s) = aux';
    zminus(:,s) = (M(:,s)' * aux)';
end
%toc

% tic
%     for c=1:nc
%         aux = 0;
%         for s=1:ns
%             aux = aux + CT(:,s) * M(c,s);
%         end
%         sumCT2dsM{c} = aux;
%     end
%
%     for s=1:ns
%         zminus{s} = zeros(m,n);
%         for c=1:nc
%             eminus{s,c} = YT(:,:,c) - sumCT2dsM{c} + CT{s} * M(c,s);
%             zminus{s} = zminus{s} + M(c,s) * eminus{s,c};
%         end
%     end
% toc

%
% if  any(abs(cell2mat(auxeminus)-cell2mat(eminus))>0.1)
%         disp('ERRORRRRRR eminus');
% end
% if  any(abs(cell2mat(auxzminus)-cell2mat(zminus))>0.1)
%         disp('ERRORRRRRR zminus');
% end


end