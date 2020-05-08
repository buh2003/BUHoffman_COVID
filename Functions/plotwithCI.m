function ax=plotwithCI(time,datay,matrixCI,real_data_length,h,varargin)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
global current_count

alpha = 0.2;

if isempty(h.Children)
    current_count = 0;
end
current_count = current_count + 1;

if current_count == 1
    acolor = [0 0.4470 0.7410];
elseif current_count == 2
    acolor = [0.8500 0.3250 0.0980];
elseif current_count == 3
    acolor = [0.9290 0.6940 0.1250];
elseif current_count == 4
    acolor = [0.4940 0.1840 0.5560];
elseif current_count == 5
    acolor = [0.4660 0.6740 0.1880];
elseif current_count == 6
    acolor = [0.3010 0.7450 0.9330];
end

CIs = [];
if size(matrixCI,2) > 1
    for i = 1:size(matrixCI,1)
        templower = 0;
        lower = 0;
        tempupper = 0;
        upper = 0;
        for j = 1:size(matrixCI,2)
            templower = matrixCI{i,j}(1);
            lower = lower+templower;
            tempupper = matrixCI{i,j}(2);
            upper = upper+tempupper;
        end
        CIs{i,1}=[lower,upper];
    end
else
    CIs = matrixCI;
end

yCIabove = [];
for i = 1:size(datay,2)
    for j = 1:size(matrixCI,1)
        yCIabove(j,i) = CIs{j,i}(2);
    end
end

yCIbelow = [];
for i = 1:size(datay,2)
    for j = 1:size(matrixCI,1)
        yCIbelow(j,i) = CIs{j,i}(1);
    end
end

if isempty(varargin) == 0
    if varargin{1}=="diff"
        yCIabove=diff(yCIabove);
        yCIbelow=diff(yCIbelow);
    end
end

fill([time(real_data_length+1:end), fliplr(time(real_data_length+1:end))],[yCIabove(real_data_length+1:end)', fliplr(yCIbelow(real_data_length+1:end)')],acolor, 'FaceAlpha', alpha,'linestyle','none');    


if ishold==0
    check=true; 
else 
    check=false;
end

hold on;ax=plot(time,datay,'Color',acolor,'linewidth',2); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end


