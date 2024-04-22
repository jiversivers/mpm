function params = changeParam(params, newValue, paramName, handles)
params.(paramName) = newValue;
params = updateParams(params);

end

