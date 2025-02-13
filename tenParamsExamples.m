rng('default');
% bfm_dir = 'C:\Users\bexlab\Dropbox\Kerri_Walter\BasalFace';
bfm_dir = 'C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace';
cd(bfm_dir); % switch to Basle Face Model Directory
[model msz] = load_model(); % load BFM
rp = defrp; % rp = render parameters for BFM
rp.phi = 0; % render front-facing view

centerPCValMu=0; % TODO - this controls the mid point of the discrimination vector
centerPCValSD=0.5; % TODO - this controls the standard deviation of the mid point of the discrimination vector

minTestLevel=findDPrimeLevel(threshEsts(pcNum),slopeEsts(pcNum),0.1,[0.1 5]); % find level with d' closest to near zero - low bound minimum disparity (0.5 pixels), high bound maximum disparity (0.5 sigma)
maxTestLevel=min([100 findDPrimeLevel(threshEsts(pcNum),slopeEsts(pcNum),4.5,[0.5 5])]); % find level with d' closest to max(5), cap at 100% contrast

centerPCVal=centerPCValMu+centerPCValSD*randn(1, nCells); % generate new rndom level of center structur values

% shape_coords = normrnd(0, 1, [199,1]); % create a random vector for the 199-D structural PCA coordinates
% tex_coords = normrnd(0, 1, [199,1]); % ...and the same for the 199-D texture PCA coordinates

for pcNum = 1:10
   
    shape_coords_control = zeros(199,1); % 0 face
    shape_coords_up = zeros(199,1); % 0 face
    shape_coords_down = zeros(199,1); % 0 face
    tex_coords = zeros(199,1); % 0 face

    rp.width=400; % rp = render parameters for BFM
    rp.height=400;
    
    %control face
    mesh1  = coef2object(shape_coords_control, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
    tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
    h = figure(1);
    display_face(mesh1, tex1, model.tl, rp);
    set(gcf, 'Color', [ 0.5 0.5 0.5 ]);
    f1 = getframe;
    img1(pcNum) = {f1.cdata};

    % render and draw first face
    shape_coords_up(pcNum) = 5; %up 5 from control

    mesh1  = coef2object(shape_coords_up, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
    tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
    h = figure(1);
    display_face(mesh1, tex1, model.tl, rp);
    set(gcf, 'Color', [ 0.5 0.5 0.5 ]);
    f1 = getframe;
    img2(pcNum) = {f1.cdata};

    % render and draw second face
    shape_coords_down(pcNum) = -5; %down 5 from control

    mesh1  = coef2object(shape_coords_down, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
    tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
    h = figure(1);
    display_face(mesh1, tex1, model.tl, rp)
    set(gcf, 'Color', [ 0.5 0.5 0.5 ])
    f1 = getframe;
    img3(pcNum) = {f1.cdata};    
end

for pcNum=1:10
    figure()
    subplot(1,3,1)
    imagesc(img1{pcNum});set(gca,'YTick',[]);set(gca,'XTick',[]);title('0');ylabel(pcNum,'FontSize',12,'FontWeight','bold');
    y = get(gca,'YLabel');set(y,'rotation',0,'VerticalAlignment','middle');
    subplot(1,3,2)
    imagesc(img2{pcNum});axis off;title('+5');
    subplot(1,3,3)
    imagesc(img3{pcNum});axis off;title('-5');
    set(gcf, 'Position',  [400, 700, 1500, 400])
end

%% Save individual faces on white background

for pcNum = 1:10
shape_coords=zeros(199,1)
shape_coords(pcNum)=5
mesh1  = coef2object(shape_coords, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
h = figure(1);
display_face(mesh1, tex1, model.tl, rp);
set(gcf, 'Color', [1 1 1 ]);
saveas(gcf,sprintf('pos%d',pcNum),'png')
end