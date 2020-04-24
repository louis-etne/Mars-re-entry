classdef FlowRegime
    enumeration
      Supersonic,
      Transonic,
      Subsonic
    end
   
    methods
        function regime = getRegime(M)
            if M >= 1.3
                regime = MFlowRegime.Supersonic;
            elseif M >= 0.8 && M < 1.3
                regime = MFlowRegime.Transonic;
            else
                regime = MFlowRegime.Subsonic;
            end
        end
    end
end