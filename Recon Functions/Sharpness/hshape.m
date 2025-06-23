function[hgram] = hshape(p)

len = p.h_ln;
gw  = p.h_sig;
dc  = p.h_dc;

x=0:len-1;
hgram=exp(-(x).^2/gw^2)+dc;
hgram=hgram/max(hgram);