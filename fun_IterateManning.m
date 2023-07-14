function [v,T,R] = fun_IterateManning(Q,B,v,M,I,maxiter,maxdiff)
            k = 0;
            v_alt = v;
                if v>0; A=Q/v; else; A=0; end
                T = A/B;
                U = 2*T+B;
                if U>0; R=A/U; else; R=0; end
                if R>0; v=M*R^.66667*I^.5; else; v=0; end
                
            while k<maxiter && abs(v_alt-v)>maxdiff
                k=k+1;
                v_alt = v;
                if v>0; A=Q/v; else; A=0; end
                T = A/B;
                U = 2*T+B;
                if U>0; R=A/U; else; R=0; end
                if R>0; v=M*R^.66667*I^.5; else; v=0; end
            end
            
        end
        