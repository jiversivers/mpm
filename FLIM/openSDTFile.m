function [out] = openSDTFile(SDTFilename)
% [out] = openSDTFile(SDTFilename)
% Load FLIM decay data, as well as associated metadata, from a becker-hickl
% *.sdt file.
%
% -------------------------------------------------------------------------
% Created by:
% Alan Woessner
% University of Arkansas
% Quantitative Tissue and Diagnostics Laboratory (www.quinnlab.org)
%
% Last edited on/by: AW, 4/21/2020
%                    AW, 5/26/2020
%                    AW, 7/1/2020
% Edit notes:
% 4/21/2020
%   -Initial Release
%       -Note: This function has only been tested with single channel,
%       single page, 512 by 512 image data. Use at your own risk. 
%
% 5/26/2020
%   -Increased compatibility to work with any sized square images. For some
%   reason, the file is shorter than expected for 1024x1024, so a band-aid
%   was used, but it works until it doesn't.
%
% 7/1/2020
%   -Removed band-aid and fixed compatibility to work with (supposedly) any size image.
%   -Consolidated outputs into a single structure variable
%   -Changed approach on reading sdt file. Makes more sense now.
% -------------------------------------------------------------------------
% Inputs:
%   SDTFilename
%       The full filename (with extension) of the .sdt file wanting to
%       load.
% -------------------------------------------------------------------------
% Outputs:
%   Structure containing:
%
%   headerInfo
%       The header portion of the .sdt file. This contains infomation about
%       when the FLIM data was acquired, as well as hex addesses associated
%       with additional portions of the file.
%
%   measurementInfo
%       The metadata associated with the hardware used to acquire the FLIM
%       data. For more information, please refer to the bh handbook for the
%       TCSPC file structure.
%
%   decayInfo
%       Portion of the file that contains hex addresses associated with the
%       decay data.
%
%   decay
%       Fluorescence decay, reported in photon counts.
% -------------------------------------------------------------------------

fID = fopen(SDTFilename);

%-- SDT Header --%
out.headerInfo.revision = fread(fID,1,'*short');
out.headerInfo.info_offset = fread(fID,1,'*long');
out.headerInfo.info_length = fread(fID,1,'*short');
out.headerInfo.setup_offs = fread(fID,1,'*long');
out.headerInfo.setup_length = fread(fID,1,'*short');
out.headerInfo.data_block_offset = fread(fID,1,'*long');
out.headerInfo.no_of_data_blocks = fread(fID,1,'*short');
out.headerInfo.data_block_length = fread(fID,1,'*ulong');
out.headerInfo.meas_desc_block_offset = fread(fID,1,'*long');
out.headerInfo.no_of_meas_desc_blocks = fread(fID,1,'*short');
out.headerInfo.meas_desc_block_length = fread(fID,1,'*short');
out.headerInfo.header_valid = fread(fID,1,'*ushort');
out.headerInfo.reserved1 = fread(fID,1,'*ulong');
out.headerInfo.reserved2 = fread(fID,1,'*ushort');
out.headerInfo.chksum = fread(fID,1,'*ushort');

%-- Information --%
out.fileInformation = reshape(fread(fID,out.headerInfo.info_length,'*char'),1,[]);

%-- Setup --%
out.setup = reshape(fread(fID,out.headerInfo.setup_length,'*char'),1,[]);

%-- Measurement Information --%
out.measurementInfo.time = fread(fID,9,'*char');
out.measurementInfo.date = fread(fID,11,'*char');
out.measurementInfo.mod_ser_no = fread(fID,16,'*char');
out.measurementInfo.measurementInfo_mode = fread(fID,1,'*short');
out.measurementInfo.cfd_ll = fread(fID,1,'*float');
out.measurementInfo.cfd_lh = fread(fID,1,'*float');
out.measurementInfo.cfd_zc = fread(fID,1,'*float');
out.measurementInfo.cfd_hf = fread(fID,1,'*float');
out.measurementInfo.syn_zc = fread(fID,1,'*float');
out.measurementInfo.syn_fd = fread(fID,1,'*short');
out.measurementInfo.syn_hf = fread(fID,1,'*float');
out.measurementInfo.tac_r = fread(fID,1,'*float');
out.measurementInfo.tac_g = fread(fID,1,'*short');
out.measurementInfo.tac_of = fread(fID,1,'*float');
out.measurementInfo.tac_ll = fread(fID,1,'*float');
out.measurementInfo.tac_lh = fread(fID,1,'*float');
out.measurementInfo.adc_re = fread(fID,1,'*short');
out.measurementInfo.eal_de = fread(fID,1,'*short');
out.measurementInfo.ncx = fread(fID,1,'*short');
out.measurementInfo.ncy = fread(fID,1,'*short');
out.measurementInfo.page = fread(fID,1,'*ushort');
out.measurementInfo.col_t = fread(fID,1,'*float');
out.measurementInfo.rep_t = fread(fID,1,'*float');
out.measurementInfo.stopt = fread(fID,1,'*short');
out.measurementInfo.overfl = fread(fID,1,'unsigned char=>char');
out.measurementInfo.use_motor = fread(fID,1,'*short');
out.measurementInfo.steps = fread(fID,1,'*ushort');
out.measurementInfo.offset = fread(fID,1,'*float');
out.measurementInfo.dither = fread(fID,1,'*short');
out.measurementInfo.incr = fread(fID,1,'*short');
out.measurementInfo.mem_bank = fread(fID,1,'*short');
out.measurementInfo.mod_type = fread(fID,16,'unsigned char=>char');
out.measurementInfo.syn_th = fread(fID,1,'*float');
out.measurementInfo.dead_time_comp = fread(fID,1,'*short');
out.measurementInfo.polarity_l = fread(fID,1,'*short');
out.measurementInfo.polarity_f = fread(fID,1,'*short');
out.measurementInfo.polarity_p = fread(fID,1,'*short');
out.measurementInfo.linediv = fread(fID,1,'*short');
out.measurementInfo.accumulate = fread(fID,1,'*short');
out.measurementInfo.flbck_y = fread(fID,1,'*int');
out.measurementInfo.flbck_x = fread(fID,1,'*int');
out.measurementInfo.bord_u = fread(fID,1,'*int');
out.measurementInfo.bord_l = fread(fID,1,'*int');
out.measurementInfo.pix_time = fread(fID,1,'*float');
out.measurementInfo.pix_clk = fread(fID,1,'*short');
out.measurementInfo.trigger = fread(fID,1,'*short');
out.measurementInfo.scan_x = fread(fID,1,'*int');
out.measurementInfo.scan_y = fread(fID,1,'*int');
out.measurementInfo.scan_rx = fread(fID,1,'*int');
out.measurementInfo.scan_ry = fread(fID,1,'*int');
out.measurementInfo.fifo_typ = fread(fID,1,'*short');
out.measurementInfo.epx_div = fread(fID,1,'*int');
out.measurementInfo.mod_type_code = fread(fID,1,'*ushort');
out.measurementInfo.mod_fpga_ver = fread(fID,1,'*ushort');
out.measurementInfo.overflow_corr_factor = fread(fID,1,'*float');
out.measurementInfo.adc_zoom = fread(fID,1,'*int');
out.measurementInfo.cycles = fread(fID,1,'*int');
out.measurementInfo.StopInfo.status = fread(fID,1,'*ushort');
out.measurementInfo.StopInfo.flags = fread(fID,1,'*ushort');
out.measurementInfo.StopInfo.stop_time = fread(fID,1,'*float');
out.measurementInfo.StopInfo.cur_step = fread(fID,1,'*int');
out.measurementInfo.StopInfo.cur_cycle = fread(fID,1,'*int');
out.measurementInfo.StopInfo.cur_page = fread(fID,1,'*int');
out.measurementInfo.StopInfo.min_sync_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.min_cfd_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.min_tac_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.min_adc_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.max_sync_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.max_cfd_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.max_tac_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.max_adc_rate = fread(fID,1,'*float');
out.measurementInfo.StopInfo.reserved1 = fread(fID,1,'*int');
out.measurementInfo.StopInfo.reserved2 = fread(fID,1,'*float');
out.measurementInfo.FCSInfo.chan = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.fcs_decay_calc = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.mt_resol = fread(fID,1,'*uint');
out.measurementInfo.FCSInfo.cortime = fread(fID,1,'*float');
out.measurementInfo.FCSInfo.calc_photons = fread(fID,1,'*uint');
out.measurementInfo.FCSInfo.fcs_points = fread(fID,1,'*int');
out.measurementInfo.FCSInfo.end_time = fread(fID,1,'*float');
out.measurementInfo.FCSInfo.overruns = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.fcs_type = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.cross_chan = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.mod = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.cross_mod = fread(fID,1,'*ushort');
out.measurementInfo.FCSInfo.cross_mt_resol = fread(fID,1,'*uint');
out.measurementInfo.image_x = fread(fID,1,'*int');
out.measurementInfo.image_y = fread(fID,1,'*int');
out.measurementInfo.image_rx = fread(fID,1,'*int');
out.measurementInfo.image_ry = fread(fID,1,'*int');
out.measurementInfo.xy_gain = fread(fID,1,'*short');
out.measurementInfo.dig_flags = fread(fID,1,'*short');
out.measurementInfo.adc_de = fread(fID,1,'*short');
out.measurementInfo.det_type = fread(fID,1,'*short');
out.measurementInfo.x_axis = fread(fID,1,'*short');
out.measurementInfo.measurementInfoHISTInfo.fida_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfo.filda_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfo.fida_points = fread(fID,1,'*int');
out.measurementInfo.measurementInfoHISTInfo.filda_points = fread(fID,1,'*int');
out.measurementInfo.measurementInfoHISTInfo.mcs_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfo.mcs_points = fread(fID,1,'*int');
out.measurementInfo.measurementInfoHISTInfo.cross_calc_phot = fread(fID,1,'*uint');
out.measurementInfo.measurementInfoHISTInfo.mcsta_points = fread(fID,1,'*ushort');
out.measurementInfo.measurementInfoHISTInfo.mcsta_flags = fread(fID,1,'*ushort');
out.measurementInfo.measurementInfoHISTInfo.mcsta_tpp = fread(fID,1,'*uint');
out.measurementInfo.measurementInfoHISTInfo.calc_markers = fread(fID,1,'*uint');
out.measurementInfo.measurementInfoHISTInfo.fcs_calc_phot = fread(fID,1,'*uint');
out.measurementInfo.measurementInfoHISTInfo.reserved3 = fread(fID,1,'*uint');
out.measurementInfo.measurementInfoHISTInfoExt.first_frame_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfoExt.frame_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfoExt.line_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfoExt.pixel_time = fread(fID,1,'*float');
out.measurementInfo.measurementInfoHISTInfoExt.scan_type = fread(fID,1,'*short');
out.measurementInfo.measurementInfoHISTInfoExt.skip_2nd_line_clk = fread(fID,1,'*short');
out.measurementInfo.measurementInfoHISTInfoExt.right_border = fread(fID,1,'*uint');
out.measurementInfo.measurementInfoHISTInfoExt.info = fread(fID,40,'unsigned char=>char');
out.measurementInfo.sync_delay = fread(fID,1,'*float');
out.measurementInfo.sdel_ser_no = fread(fID,1,'*ushort');
out.measurementInfo.sdel_input = fread(fID,1,'unsigned char=>char');
out.measurementInfo.mosaic_ctrl = fread(fID,1,'unsigned char=>char');
out.measurementInfo.mosaic_x = fread(fID,1,'unsigned char=>char');
out.measurementInfo.mosaic_y = fread(fID,1,'unsigned char=>char');
out.measurementInfo.frames_per_el = fread(fID,1,'*short');
out.measurementInfo.chan_per_el = fread(fID,1,'*short');
out.measurementInfo.mosaic_cycles_done = fread(fID,1,'*int');
out.measurementInfo.mla_ser_no = fread(fID,1,'*ushort');
out.measurementInfo.DCC_in_use = fread(fID,1,'*uchar');
out.measurementInfo.dcc_ser_no = fread(fID,12,'unsigned char=>char');
out.measurementInfo.TiSaLas_status = fread(fID,1,'*ushort');
out.measurementInfo.TiSaLas_wav = fread(fID,1,'*ushort');
out.measurementInfo.AOM_status = fread(fID,1,'unsigned char=>char');
out.measurementInfo.AOM_power = fread(fID,1,'unsigned char=>char');
out.measurementInfo.ddg_ser_no = fread(fID,8,'unsigned char=>char');
out.measurementInfo.prior_ser_no = fread(fID,1,'*int');
out.measurementInfo.mosaic_x_hi = fread(fID,1,'unsigned char=>char');
out.measurementInfo.mosaic_y_hi = fread(fID,1,'unsigned char=>char');
out.measurementInfo.reserve = fread(fID,12,'unsigned char=>char');

%-- Data block header --%

out.decayInfo.data_offs_ext = fread(fID,1,'*uchar');
out.decayInfo.next_block_offs_ext = fread(fID,1,'*uchar');
out.decayInfo.data_offs = fread(fID,1,'*ulong');
out.decayInfo.next_block_offs = fread(fID,1,'*ulong');
out.decayInfo.block_type = fread(fID,1,'*ushort');
out.decayInfo.meas_desc_block_no = fread(fID,1,'*short');
out.decayInfo.lblock_no = fread(fID,1,'*ulong');
out.decayInfo.block_length = fread(fID,1,'*ulong');

%-- Decay data --%
decay = fread(fID,'ushort=>uint16');
fclose(fID);

blankDecay = zeros(out.measurementInfo.image_y,out.measurementInfo.image_x,out.measurementInfo.adc_re,'uint16');
counter = 0;
for i = 1:length(decay)/numel(blankDecay)
    tmpDecay = blankDecay;
    for j = 1:double(out.measurementInfo.adc_re)-1
        tmpDecay(:,:,j) = rot90(reshape(decay(counter+j:double(out.measurementInfo.adc_re):counter+j+numel(blankDecay)-1),out.measurementInfo.image_x,out.measurementInfo.image_y),1);
    end
out.decay.(strcat('Channel_',int2str(i))) = tmpDecay;
counter = counter+numel(blankDecay); 
end

end