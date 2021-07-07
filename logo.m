format shortG

%% Load the image
I = imread('logo.png');
I= im2bw(I);
dim = size(I);
boundaries = bwboundaries(I);
numboundaries = size(boundaries)
figure
imshow(I);
hold on;
for k=1:numboundaries(1)
    b= boundaries{k};
    plot(b(:,2),b(:,1),'g','Linewidth',3);
end
hold off;
shield = boundaries{6};
scimatar1 = boundaries{2};
scimatar2 = boundaries{3};
scimatar3 = boundaries{4};
crown = boundaries{5};
wholelogo = [shield;scimatar1;scimatar2;scimatar3;crown];

    %% Choose an element:
    element = crown;
    complexelement = element(:,2) + i*element(:,1);

    %% Discrete Fourier Transform: n.b. normalise here.
    N=size(complexelement,1);
    coefficients = 1/N * exp(-i*2*pi*[0:N-1]'*[0:N-1]/N) * complexelement;

    %% As frequency, circle radii and phase shifts
    circles = [[0:N-1]' abs(coefficients) atan2(imag(coefficients),real(coefficients))];

    %% Arrange by most significant mode
    % These are the values we need to feed into the reconstruction.
    % These should be adjusted to so that the frequencies are as small as
    % possible (allowing negative frequencies) by subtracting N, otherwise
    % the interpolation is unnecessarily complicated!
    % Todo: automate this!
    sortedcircles = sortrows(circles,2,'descend');

    %csvwrite('output.csv',sortedcircles([1:8],[2,1,3]))

%% Reconstitute function
t = [0:N-1];
interp = coefficients'*exp(i*2*pi*[0:N-1]'*t/N);
plot(interp)


%% Trim some modes:
%% Number of modes: k
k = 9;
interp = coefficients(1:k)'*exp(i*2*pi*[0:k-1]'*t/N);
plot(interp)

%% Animate!
F(N) = struct('cdata',[],'colormap',[]);
for k = 1:N
    interp = coefficients(1:k)'*exp(i*2*pi*[0:k-1]'*t/N);
    plot(interp)
    drawnow
    F(k) = getframe;
end
movie(F)

%% Sort the modes first:
sortedcoefficients = [[0:N-1]' abs(coefficients) coefficients];
sortedcoefficients = sortrows(sortedcoefficients,2,'descend');

%% Trim some modes:
%% Number of modes: k
k =3;
interp = sortedcoefficients(1:k,3)'*exp(i*2*pi*sortedcoefficients(1:k,1)*t/N);
plot(interp)

%% DFT Interpolation
t = 0:0.001:1;
interp = sortedcoefficients(1:k,3)'*exp(i*2*pi*sortedcoefficients(1:k,1)*t);
plot(interp)


%% Animate!
F(N) = struct('cdata',[],'colormap',[]);
for k = 1:N
    interp = sortedcoefficients(1:k,3)'*exp(i*2*pi*sortedcoefficients(1:k,1)*t/N);
    plot(interp)
    drawnow
    F(k) = getframe;
end
movie(F)
