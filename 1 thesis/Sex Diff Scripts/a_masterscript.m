%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%                           MASTER SCRIPT                            %
%                                 -                                  %
%                  Thesis: Sex Differences in Neglect                %
%                            Tamara Keßler                           %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This master script ONLY includes custom scripts, but none for:
% 1: Conversion DICOM to NII (-> MRIConvert.exe)
% 2: Lesion Delineation (-> Clusterize Toolbox)
% 3: Image Pre-Processing, incl. Coregistration & Normalisation
%    (-> SPM12 & Clinical Toolbox)
% 4: Lesion Quantification Toolkit

%% Setting up

clear; clc;

% Set up struct for the relevant directories
dir = struct;

% Home directory for the scripts
dir.home = "C:\Users\Sperber\Desktop\Sex Diff Scripts";
cd(dir.home)

% Directories for the data
dir.scans = "D:\Tamara\Scans";
dir.lesionmaps = "D:\Tamara\LesionMaps";
dir.disconmaps = "D:\Tamara\DisconnectionMaps";

%% (Fix fragmented scans)

% If patient's scans have been fragmented for some reason
% - meaning their full brain scan is split into 2-3 different
% NII files: Use this script to combine them into a single NII

concat_scans

%% (Add hypo- & hyperdense lesion maps)

% If patient suffered an infarct, as well as a haemorrhage and
% now exhibits hypo- & hyperdense regions: Use this script to
% combine them into a single lesion map, while controlling for overlaps

add_lesion_maps

%% Check Bounding Box

% Use this script to check if all scans/lesion maps have the same
% bounding box (181x217x181)

Check_Bounding_Box_and_Max_Value_new

%% Calculate lesion size

% Use this script to calculate the lesion size for every patient

get_lesion_size

%% Fix Origin

% Use this script to change the origin of every patient's lesion map
% from [-90, -126, -72] to [-90, -125, -71]

FixOrigin_nii_edited

%% Mask Ventricles

% Use this script to mask (aka cut out) the ventricles and cerebellum
% from the lesion maps, to make them fit the ch2bet template better

Mask_Ventricles_edited

%% Prep for NiiStat

convenience
rename_files_4niistat

%% [NiiStat]

%% Move Origin for LQT

% Use this toolbox to change the lesion maps' dimensions from (181x217x181)
% to (182x218x182), which is required for the Lesion Quantification Toolbox

move_nii_origin_edited


%% [LQT]

%% Check LQT results before further steps

sanity_check

%% Create Visualisations

% (Create mat files of the Disconnection Maps using nii_nii2mat)
move_matfiles
create_overlays_and_average_disconnections

%% Prep for ROI-to-ROI GLM

move_disconnectionmaps
adjust_pnums

%% ROI-to-ROI GLM

massunivariate_disc_GLM_edited_morePs
eval_parcel_disconns
