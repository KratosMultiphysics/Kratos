% clear all;
% close all;


a = importdata('_list_demBig_1.post.lst');

Variable1 = '"DISPLACEMENT"'
Variable_Solution1 = {};
solution = [0,0.0,0.0,0.0];
Variable_Solution1{end+1} = solution;

Variable2 = '"VELOCITY"'
Variable_Solution2 = {};
solution = [0,0.0,0.0,0.0];
Variable_Solution2{end+1} = solution;

Variable3 = '"ANGULAR_VELOCITY"'
Variable_Solution3 = {};
solution = [0,0.0,0.0,0.0];
Variable_Solution3{end+1} = solution;

for i = 2:length(a)
   filename = string(a(i));
   filename = filename{1};
   fileformat = filename(end-2:end)
   
   if (fileformat == 'res')
       file_line = fopen(filename)
       
       tline = fgetl(file_line);
       tline = fgetl(file_line);
       while ischar(tline)
           b = strsplit(tline);
           if strcmp(b{1},'Result')
               solution_step = str2double(b{4});
               if (b{5} == 'Vector')
                   tline = fgetl(file_line);
                   tline = fgetl(file_line);
                   while ~strcmp(tline,'End Values')
                       solution_values = strsplit(tline)
                       id = str2num(solution_values{1});
                       solution = [solution_step,
                           str2double(solution_values{2}),
                           str2double(solution_values{3}),
                           str2double(solution_values{4})];
                       
                       if strcmp(b{2}, Variable1)
                           Variable_Solution1{end+1} = solution;
                       end
                       if strcmp(b{2}, Variable2)
                           Variable_Solution2{end+1} = solution;
                       end
                       if strcmp(b{2}, Variable3)
                           Variable_Solution3{end+1} = solution;
                       end
                       tline = fgetl(file_line);
                   end
               end
           end
           tline = fgetl(file_line);
       end
       fclose(file_line);
   elseif(fileformat == 'msh')
       
   end
end


X1 = [];
Y1 = [];
Y2_analytical = [];
Error1 = [];
for i = 1:length(Variable_Solution1)
    X1(end+1) = Variable_Solution1{i}(1);
    Y1(end+1) = Variable_Solution1{i}(2);
    Y2_analytical(end+1) = velocity(Variable_Solution1{i}(2));
%     Error1(end+1) = (Y1(end)-Y1_analytical(end))/Y1_analytical(end);
end

X2 = [];
Y2 = [];
% Y2_analytical = [];
% Error2 = [];
for i = 1:length(Variable_Solution2)
    X2(end+1) = Variable_Solution2{i}(1);
    Y2(end+1) = Variable_Solution2{i}(2);
%     Y2_analytical(end+1) = velocity(Variable_Solution2{i}(1));
%     Error2(end+1) = (Y2(end)-Y2_analytical(end))/Y2_analytical(end);
end

% X3 = [];
% Y3 = [];
% Y3_analytical = [];
% Error3 = [];
% for i = 1:length(Variable_Solution3)
%     X3(end+1) = Variable_Solution3{i}(1);
%     Y3(end+1) = Variable_Solution3{i}(3);
%     Y3_analytical(end+1) = angular_velocity(Variable_Solution3{i}(1));
%     Error3(end+1) = (Y3(end)-Y3_analytical(end))/Y3_analytical(end);
% end

% figure;
hold on;
% xlim([0 0.25])
% ylim([0 13])
xlabel('time [s]') % x-axis label
ylabel('displacement [m], velocity [m/s], angular velocity [rad/s]') % y-axis label
% plot(X1, Y1_analytical, 'color', [0,101,189]/255)
plot(X1, Y1, 'color', [0,101,189]/255)
plot(X1, Y2_analytical, 'color', [88,88,90]/255)
plot(X2, Y2, 'color', [0.0,82,147]/255)
% plot(X3, Y3_analytical, 'color', [0.0,51,89]/255)
% plot(X3, Y3, 'color', [88,88,90]/255)
legend('displacements analytical','displacements','velocity analytical','velocity','angular velocity analytical','angular velocity','Location','east')

% figure;
% hold on;
% xlim([0 1])
% ylim([-5e-3 5e-3]) 
% xlabel('time [s]') % x-axis label
% ylabel('error [%]') % y-axis label
% plot(X1,Error1, 'color', [217,218,219]/255)
% plot(X2,Error2, 'color', [156,157,159]/255)
% plot(X3,Error3, 'color', [88,88,90]/255)
% legend('displacements','velocity','angular velocity','Location','east')

% function av = angular_velocity(time)
% R = 0.3;
% mu = 0.3;
% m = 100*(4/3)*pi*R*R*R;
% g = 9.81;
% v_0 = 5;
% tc = (2*v_0)/(7 * mu* g);
% if(time <= tc)
%  av = (R*mu*m*g*time)/((2 * m *R*R)/5);
% else
%  av = (5 * v_0)/(7*R);
% end
% end

function v = velocity(distance)
r = 0.05172;
% mu = 0.3;
% m = 100*(4/3)*pi*R*R*R;
g = 9.81;
% v_0 = 5;
% tc = (2*v_0)/(7 * mu* g);
distance
 v = sqrt(2*g*r*(1-cos(distance)));
v
end

function x = displacement(time)
R = 0.3;
mu = 0.3;
m = 100*(4/3)*pi*R*R*R;
g = 9.81;
v_0 = 5;
tc = (2*v_0)/(7 * mu* g);
x = time / sqrt(0.05172/g)
% if(time <= tc)
%  x = v_0 * time - mu * g * time * time / 2;
% else
%  x = (12 * v_0 *v_0)/(49*mu*g) + 5*v_0*(time-tc)/7;
% end
end