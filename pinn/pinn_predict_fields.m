function pred = pinn_predict_fields(net, xyt)
% pinn_predict_fields
%
% Evaluates the PINN at input points xyt = [x, y, t].
%
% INPUTS:
%   net - trained or untrained dlnetwork
%   xyt - numeric array of size [N, 3]
%
% OUTPUT:
%   pred - structure containing named predicted fields

xyt_dl = dlarray(xyt.', 'CB');
out = predict(net, xyt_dl);
out = extractdata(out).';

pred = struct();
pred.rho = out(:,1);
pred.ux  = out(:,2);
pred.uy  = out(:,3);
pred.uz  = out(:,4);
pred.p   = out(:,5);
pred.Bx  = out(:,6);
pred.By  = out(:,7);
pred.Bz  = out(:,8);

end
