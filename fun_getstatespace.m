function [m_A,m_B,m_C,m_D] = fun_getstatespace(n,v_kpu,v_b1,v_b2,v_a1,v_a2)

% This function computes the state space representation matrices A, B, C, D
% given the canonical transfer function form
%
% Each unit is represented by a generic 2nd order transfer
% function: P(s) = kpu*(1 + b1*s + b2*s^2)/(1 + a1*s + a2*s^2)
%
% A = [0 1; -1/a2 -a1/a2]
% B = [0; 1/a2]
% C = kpu*[(1 - b2/a2) (b1 - b2*a1/a2)]
% D = kpu*b2/a2

m_A = zeros(2*n,2*n);
m_B = zeros(2*n,1);
m_C = zeros(n,2*n);
m_D = zeros(n,1);

for i = n:-1:1
   if v_a2(i)>0
      m_A(2*(i-1)+1:2*i,2*(i-1)+1:2*i) = [0 1; -1/v_a2(i) -v_a1(i)/v_a2(i)]; 
      m_B(2*i,:) = 1/v_a2(i);
      m_C(i,2*(i-1)+1:2*i) = [(-v_kpu(i)-v_kpu(i)*v_b2(i)/v_a2(i)) (-v_kpu(i)*v_b1(i)-v_kpu(i)*v_b2(i)*v_a1(i)/v_a2(i))];
      m_D(i,:) = v_kpu(i)*v_b2(i)/v_a2(i);
   else % if a2 = 0 (1st order transfer function)
      m_A(2*(i-1)+1:2*i,2*(i-1)+1:2*i) = [0 0; 0 -1/v_a1(i)];  
      m_B(2*i,:) = 1/v_a1(i);
      m_C(i,2*(i-1)+1:2*i) = [0 (-v_kpu(i)+v_b1(i)*v_kpu(i)/v_a1(i))];
      m_D(i,:) = -v_kpu(i)*v_b1(i)/v_a1(i);
   end
end

