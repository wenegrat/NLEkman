qvec = 0:0.01:6;
l = -freq/2:freq/cuts:freq/2;

fl = .5*cos(2*pi.*l);

flfft = fft(fl);
% qvec = abs(flfft);
%%

for q = qvec;
    
    n=10;
    cuts = 500;
    freq = 2;
    kn2=[];
    for jcut=-freq/2:freq/cuts:freq/2
        A= zeros(2*n+1);
        for j=1:2*n+1
            A(j,j) = (n+1-j+jcut).^2;
        end
        for j=1:2*n+1-freq
            A(j+freq, j) = q;
            A(j,j+freq) = q;
        end
        [V, W] = eig(A);
        kn2=[kn2; diag(W)];
    end
%     plot((real(kn2)./sqrt(abs(real(kn2)))), (imag(kn2)+sqrt(q)), 'k.'); hold on; axis([-20 20 0 6])
     plot(real(kn2), imag(kn2)+sqrt(q), 'k.'); hold on; axis([-20 20 0 6])

end  

xlabel('k_n^2'); ylabel('q');
hold off