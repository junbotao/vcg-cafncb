classdef MicroCounter < handle
    properties
        succ = 0;
        att  = 0;
    end
    methods
        function reset(obj)
            obj.succ = 0;
            obj.att  = 0;
        end
    end
end
