function [t_out,y_out] = fastode45_MEE(tspan,y0,c, Thr,rho,SI2CAN) %#codegen
global atol rtol pow A B E

t0 = tspan(1);
tfinal = tspan(2);
f0 = StCO_Dynamics_MEE(t0,y0, c, Thr,rho,SI2CAN);
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = 0.1*(tfinal-t0);
tnew = 0;
ynew = y0;
t = t0;
y = y0;
Neq = 14;                                                                   %<----------------
f = zeros(Neq,7,'double'); % the first number is number of equations
hmin = 16*eps(t);
htspan = tfinal-t0;

% Compute an initial step size h using y'(t).
absh = min(hmax, htspan);
rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);

if absh * rh > 1
    absh = 1 / rh;
end
absh = max(absh, hmin);

f(:,1) = f0;

% THE MAIN LOOP

done = false;
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
      hA = h * A;
      hB = h * B;
      f(:,2) = StCO_Dynamics_MEE(t+hA(1),y+f*hB(:,1), c, Thr,rho,SI2CAN);
      f(:,3) = StCO_Dynamics_MEE(t+hA(2),y+f*hB(:,2), c, Thr,rho,SI2CAN);
      f(:,4) = StCO_Dynamics_MEE(t+hA(3),y+f*hB(:,3), c, Thr,rho,SI2CAN);
      f(:,5) = StCO_Dynamics_MEE(t+hA(4),y+f*hB(:,4), c, Thr,rho,SI2CAN);
      f(:,6) = StCO_Dynamics_MEE(t+hA(5),y+f*hB(:,5), c, Thr,rho,SI2CAN);
      
      tnew = t + hA(6);
      if done
          tnew = tfinal;   % Hit end point exactly.
      end
%       h = tnew - t;      % Purify h.
      
      ynew = y + f*hB(:,6);
      f(:,7) = StCO_Dynamics_MEE(tnew,ynew, c, Thr,rho,SI2CAN);
      % Estimate the error.
      err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
      
      
      % Accept the solution only if the weighted error is no more than the
      % tolerance rtol.  Estimate an h that will yield an error of rtol on
      % the next step or the next try at taking this step, as the case may be,
      % and use 0.8 of this value to avoid failures.
      if err > rtol                       % Failed step
          if absh <= hmin
              
              t_out = 0;
              y_out = zeros(Neq,1);
%             sprintf('tolerance not satisfied');
              return;
          end
          
          if nofailed
              nofailed = false;
              absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
              
          else
              absh = max(hmin, 0.5 * absh);
          end
          h = absh;
          done = false;
      else
          
          break;
      end
  end
             
  if done
    
    break;
  end
  
  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  
  f(:,1) = f(:,7);  % Already have f(tnew,ynew)
 
end
t_out = tnew;
y_out = ynew;

end



  
