function params = updateParams(params)
    params.w = 2*pi*params.f;
    if params.calc 
        params.dt = 1/(1e6*params.f*params.T);
    end
end