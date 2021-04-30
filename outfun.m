function stop = outfun(x,optimValues,state,history)
    stop = false;
    history.add(optimValues.funccount,optimValues.fval);
end