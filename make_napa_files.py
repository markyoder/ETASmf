import vc_parser

napa_sections = vc_parser.napa_region_section_filter['filter']
# and these can just mix with the existing EMC data.
A= vc_parser.get_EMC_CFF(sections=napa_sections, f_out_pattern='data/VC_CFF_timeseries_section_%d.npy')

# now, run some statistics:
b=vc_parser.waiting_time_figs(section_ids=napa_sections, m0=7.0, output_dir='VC_Napa_m7_WT_figs')

b=vc_parser.recurrence_figs(section_ids=napa_sections, m0=7.0, output_dir='VC_Napa_m7_recurrence_figs')

b=vc_parser.waiting_time_figs(section_ids=napa_sections, m0=6.0, output_dir='VC_Napa_m6_WT_figs')

b=vc_parser.recurrence_figs(section_ids=napa_sections, m0=6.0, output_dir='VC_Napa_m6_recurrence_figs')
