classdef histclass < handle
    properties 
        count= [];
        fval= [];
    end
    methods 
        function add(obj,newx,newfval)
            obj.count= newx;
            obj.fval= [obj.fval  newfval];
        end
    end
end