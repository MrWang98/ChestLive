function ewt = EWT_Single_filter(ff,mfb,l)

ewt=real(ifft(conj(mfb).*ff));
ewt=ewt(l:end-l);
