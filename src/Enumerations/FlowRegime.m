classdef FlowRegime
    enumeration
      Supersonic,
      Transonic,
      Subsonic
    end
   
    methods
        function regime = getRegime(M)
            if M >= 1.3
                regime = FlowRegime.Supersonic;
            elseif M >= 0.8 && M < 1.3
                regime = FlowRegime.Transonic;
            else
                regime = FlowRegime.Subsonic;
            end
        end
    end
end