function [ReceivedSamples]=UWOC_channel_ae(aU,bU,cU,OutputSamples,h,SNR,Nt,Nr,FlipFlag,Eb_N0_dB)

Fov = 30; % Receiver's FOV 
%constants %variables of medium 
Link_Len=30;%m
z_min = 0; z_max = Link_Len; 
% refractive index according to McNiel approximation.
lamda = 470; %nm wavelength(as per paper)
% temperature
T=13.07;
% salinity
S=33.4160;
n_waterc = 1.3247 - 2.5e-6*T^2 + S*(2e-4 - 8e-7*T) + 3300/lamda^2 - 3.2e7/lamda^4;
    
n_water = min(1.333,n_waterc); 
n_air = 1; 
n_bottom = 1.45; 
roulete_constand = 1; % This is for roulette functionality 
thetaTr = pi/4; 
phiTr = pi/3; 
positionTr = [0,0,0]; 
%variables of receiver 
thetaRe = pi/2; 
positionRe = [20,0,0]; 
radiusRe = 0.2; 
%variables of photons 
         NumSamples = size(OutputSamples,2);

N = 1e3; 
weight_thr = 10^(-4); %initialize arrays for photons 
% The total distance that each photon traveled
photons_distance = zeros(N,1);

photons_weight = zeros(N,1); 
photons_time = zeros(N,1); % The time that each photon reached the receiver 
photons_received = 0; 
total_intensity = 0; %speed in the water 
speed_in_space = 3*1e8; 
v = speed_in_space / n_water; 
countL=1;
for i=1:N %initialize photon 
    [ thetaPh, phiPh, positionPh, directionPh, weightPh ] = initialize_photon( thetaTr, phiTr, positionTr );
    reached = false; 
    total_distance = 0; 
    total_time = 0; 
%     while ~reached && weightPh > 0
% for iop=1:5
        %move the photon 
        d = - log(rand)/cU; 
        positionPh_new = move_photon(d, positionPh, directionPh);
        
        %We check if the photon exceeds the medium 
        if positionPh_new(3) > z_max 
            % If the photon gets to the air 
            % Photon's new position is on the surface 
            z_new = z_max; 
            d = (z_max - positionPh(3))/directionPh(3);%distance 
            x_new = positionPh(1) + directionPh(1)*d; 
            y_new = positionPh(2) + directionPh(2)*d; 
            positionPh_new = [x_new, y_new, z_new]; 
            t = d/v; 
            total_distance = total_distance + d; 
            total_time = total_time + t;
            %computation of Fresnel reflection coefficient 
            %refraction angle 
            thetaRefr = asin(n_water*sin(thetaPh)/n_air); 
            % Fresnel coefficient 
            fresnel_co = 1/2 * (((sin(thetaPh-thetaRefr))^2) /((sin(thetaPh+thetaRefr))^2) + ((tan(thetaPh-thetaRefr))^2) /((tan(thetaPh+thetaRefr))^2)); 
            % We check if we have internal reflection 
            rand_var = rand*1e-3; 
            if rand_var < fresnel_co 
                %internal reflection 
                directionPh(3) = -directionPh(3); % z_new = -z 
            else
                %in the air 
                weightPh = 0; 
            end
        elseif positionPh_new(3) < z_min 
            %If the photon reached the bottom 
            %Photon's new position is on the bottom 
            z_new = z_min;
            d = (z_min - positionPh(3))/directionPh(3); 
            x_new = positionPh(1) + directionPh(1)*d; 
            y_new = positionPh(2) + directionPh(2)*d; 
            positionPh_new = [x_new, y_new, z_new];
            t = d/v; 
            total_distance = total_distance + d; 
            total_time = total_time + t;
            % refraction angle 
            thetaRefr = asin(n_water*sin(thetaPh)/n_bottom); 
            % Fresnel coefficient 
            fresnel_co = 1/2 * (((sin(thetaPh-thetaRefr))^2) /((sin(thetaPh + thetaRefr))^2) + ((tan(thetaPh - thetaRefr))^2) /((tan(thetaPh + thetaRefr))^2)); 
            % At this point z = z_min. But we still need to check if the photon is received. 
            [reached, positionPh_received, out_of_fov] = received(positionRe, radiusRe, positionPh, directionPh, positionPh_new, Fov); 
            if (reached && out_of_fov==false)
                positionPh_new = positionPh_received; 
                photons_received = photons_received + 1; 
                total_intensity = total_intensity + weightPh;
                total_intensityN(countL)=weightPh;
                countL=countL+1;
                photons_time(i) = total_time; 
                break; 
                % since the photon is received, there is no need to calculate absorption and scattering 
            elseif (reached==false && out_of_fov) 
                positionPh_new = positionPh; 
            end
            %We check if we have internal reflection 
            rand_var = rand*1e-3; 
            if rand_var < fresnel_co 
                %internal reflection 
                directionPh(3) = -directionPh(3); % z_new = -z 
            else
                %in the sand -> it is lost 
                weightPh = 0;
            end
        else
            %The photon is in the water 
            % Check if the photon is received 
            [reached, positionPh_received, out_of_fov] = received(positionRe, radiusRe, positionPh, directionPh, positionPh_new, Fov); 
            distance = norm(positionPh-positionPh_received)^2;
            t = distance/v; total_distance = total_distance + distance; 
            total_time = total_time + t;
            if (reached && out_of_fov==false) 
                positionPh_new = positionPh_received;
                photons_received = photons_received + 1; 
                total_intensity = total_intensity + weightPh;
                total_intensityN(countL)=weightPh;
                countL=countL+1;
                photons_time(i) = total_time;
                break; % since the photon is received, there is no need to calculate absorption and scattering 
            elseif (reached==false && out_of_fov) 
                positionPh_new = positionPh; 
            end
            
            % The photon is in the water and not received. 
%             So the Monte Carlo continues 
            %absorption 
            weightPh_new = absorb(weightPh, aU, cU); 
            weightPh = weightPh_new; % We check the photon's weight 
            if (weightPh < weight_thr) 
                %roulette 
                propability_of_survival = 1 / roulete_constand; 
                x = rand; 
                if (x <= propability_of_survival) 
                    weightPh = roulete_constand * weightPh; 
                else
                    weightPh = 0; 
                    break; 
                end
            end % scattering
%             scattering method
          ScaterMethod='TTHG';%HG
           g = 0.924; % HG asymmetry parameter %variables of transmitter 
 [thetaPh_new, phiPh_new, directionPh_new,FadingS] = scattering(g,aU, bU,cU, directionPh,ScaterMethod); 
            thetaPh = thetaPh_new; 
            phiPh = phiPh_new; 
            directionPh = directionPh_new*FadingS; 
        end
        positionPh = positionPh_new; 
%     end %in the end, before moving to the next photons we save the followings
    photons_direction(i,:) = directionPh;
    photons_weight(i) = weightPh;
    photons_distance(i) = total_distance;
end
intensity = 10*log10(total_intensity/N); 
photons_percentage = (photons_received/N)*100; 
chanFad=aU+bU+cU;

optics.type = 'CPC';
optics.params = Fov*pi/180; % In this case corresponds to 30 degrees FOV FWHM
optics.area = pi/4*(5e-3)^2; % 5mm diameter receiver
optics.orientation = [0,0,-1]; % pointing vector
optics.n = 1.45;

% CPC's cosine of FWHM FOV/2. This is used to check if a ray contributes 
% effectively or not
cos_half_fov = cos(optics.params/2);

% Lens gain (CPC gain)
lens_gain = (optics.n/sin(optics.params/2))^2;

 t = photons_distance*n_water/3e8;
 counter = 1;
 for kl=1:N
    impact_cosine = -photons_direction(kl,:)*optics.orientation';
    
    % Entry condition check (ray is within rx's FOV)
%     if (cos_half_fov < impact_cosine)
        time(counter) = t(kl);
%         channel impulse response
        h_t(counter) = impact_cosine*optics.area*lens_gain;
         counter = counter + 1;
%     end
    
 end
%figure(1);
%subplot(211)
%stem(real(h_t));
%title('UWOC Response');     
h_t=max(h_t/max(h_t));

UWOC_receivedSig=OutputSamples*h_t;
        Nf = size(UWOC_receivedSig,1);
        % channel noise addition
        sigma = (2e-1*chanFad)+(1/sqrt(SNR)); % noise variance
        for n = 1:Nf
            UWOC_receivedSig(n,:) = UWOC_receivedSig(n,:) / sqrt(mean(abs(UWOC_receivedSig(n,:)).^2)); % normalize power
        end
        if size(UWOC_receivedSig,2)==1
            UWOC_receivedSig = UWOC_receivedSig';
            FlipFlag = 1;
        end
        
        if Nt==1 && (Nr==1 || Nr==2) %(receive diversity)
            ReceivedSamples = h.*UWOC_receivedSig + (10^(-Eb_N0_dB/20)*sigma*randn(Nr,size(UWOC_receivedSig,2)) + j*sigma*randn(Nr,size(UWOC_receivedSig,2)));
        elseif Nt==2 && Nr==1 %(transmit diversity)
            ReceivedSamples = sum(h.*UWOC_receivedSig,1) + (10^(-Eb_N0_dB/20)*sigma*randn(1,size(UWOC_receivedSig,2)) + j*sigma*randn(1,size(UWOC_receivedSig,2)));
        elseif Nt==2 && Nr==2 %MIMO
            ReceivedSamples = squeeze(sum(h.*UWOC_receivedSig,2)) + (10^(-Eb_N0_dB/20)*sigma*randn(Nr,size(UWOC_receivedSig,2)/Nt) + j*sigma*randn(Nr,size(UWOC_receivedSig,2)/Nt));
        end
        


                
        
   if FlipFlag
      ReceivedSamples = ReceivedSamples';
   end
%figure(1)
%subplot(212)
%if Nt>1 && Nr>1
%    surf(real(reshape(UWOC_receivedSig(:,:,1:50),Nt*Nr,[])))
%elseif Nt==1 && Nr==1
%    surf(real(reshape(UWOC_receivedSig(:,1:16),2,[])))
%else
%    surf(real(UWOC_receivedSig(:,1:50)))
%end
%xlabel('No of Samples');ylabel('No of MIMO antenna');zlabel('Channel Coefficient');
%title('Channel Signal');
