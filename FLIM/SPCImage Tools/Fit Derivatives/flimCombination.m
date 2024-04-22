function flimData = flimCombination(flimData)

flimData.a1a2ratio = flimData.a1./flimData.a2;
flimData.tm = (flimData.a1.*flimData.t1 + flimData.a2.*flimData.t2)./(flimData.a1+flimData.a2);

end
