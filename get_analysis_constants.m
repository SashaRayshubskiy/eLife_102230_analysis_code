function aconstants = get_analysis_constants

aconstants.TRIAL_TYPE_CNT = 3;
aconstants.LEFT = 1;
aconstants.RIGHT = 2;
aconstants.BOTH = 3;

aconstants.task_str = {'LeftOdor', 'RightOdor', 'BothOdor'};

aconstants.LEFT_CLR =  rgb('FireBrick');
aconstants.RIGHT_CLR = rgb('SeaGreen');
aconstants.BOTH_CLR = rgb('DarkBlue');

aconstants.LEFT_CLR_SINGLE =  rgb('LightSalmon');
aconstants.RIGHT_CLR_SINGLE = rgb('PaleGreen');
aconstants.BOTH_CLR_SINGLE = rgb('LightBlue');

aconstants.clr_by_type = { aconstants.LEFT_CLR, aconstants.RIGHT_CLR, aconstants.BOTH_CLR };
aconstants.single_clr_by_type = { aconstants.LEFT_CLR_SINGLE, aconstants.RIGHT_CLR_SINGLE, aconstants.BOTH_CLR_SINGLE };

aconstants.VEL_FWD = 1;
aconstants.VEL_SIDE = 2;
aconstants.VEL_YAW = 3;

aconstants.order    = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('DarkSlateGray'); rgb('Purple'); rgb('Goldenrod'); rgb('Indigo'); rgb('Olive'); rgb('Cyan'); rgb('MediumVioletRed')];

end

