import pandas as pd

df = pd.read_csv("DistanceMatrix.csv", index_col=0)

clusters = {
    1:  ['1ag9_1_B','2cnh_1_A','1f4n_2_A','1fz7_4_AE','1g8f_7_A','1jn9_1_AB','1uww_4_A','2afb_7_B','2cyy_3_A','2f1w_2_A','2gjv_5_AB','2rld_1_A','2w3o_1_B','2xn6_2_A','2zs0_4_BC','3b3d_2_B','3bs6_2_A','3fmg_1_A','3gv4_4_A','3jru_2_B','3n08_2_A','3p85_1_A','3s6f_1_A','3s82_1_A','4bcu_1_A','4f53_2_A','4fgc_3_DE','4i9f_2_B','4ii3_4_C','4jbe_4_AB','4jk4_4_B','4mbo_1_A','4p3q_1_A','5a55_5_A','5b66_74_AEFV','5buo_4_A','5cyb_2_A','5cyb_3_A','5eto_3_A','5f47_2_B','5nq2_1_A','5uo9_5_AC','5v4d_1_E','6bn2_3_A','6e54_1_A'],
    2:  ['1bob_1_A','4k3l_8_B','1v6c_4_B','2bko_5_A','2g0y_2_A','2pkt_2_A','2y5p_5_C','3b3d_4_A','3bvc_1_B','3e86_3_A','3goa_3_B','3mq6_6_CN','3qxg_6_B','3r4i_6_D','3rbz_12_B','3rmk_14_F','3t66_1_A','4ecg_4_A','4lfy_5_AB','4lfy_1_AB','4onu_2_A','4qel_2_A','4zgf_1_A','5abq_4_A','5e2f_3_A','5gmf_4_A','5hzl_3_B','5i2n_1_B','5nsa_1_A','5o7u_3_B','6eq1_3_B','6gsz_5_A','6hun_1_A','7c2g_4_G'],
    3:  ['1d0b_2_A','4hhr_9_A','1v3w_1_A','2bko_1_A','2d4c_1_A','2i52_6_C','2id3_7_B','2w4y_4_A','3b40_4_A','4hoj_2_A','4jk4_2_A','4yis_3_B'],
    4:  ['1en7_1_B','2xts_10_AD','1snn_1_A','1yax_1_A','1zmr_2_A','2fe1_1_A','2rji_1_A','7c23_3_AB','2yll_2_A','3cv1_2_A','3de8_4_B','3dzc_1_A','3e11_2_B','3fmg_2_A','3g5r_1_A','3ibv_1_A','3ljk_1_A','3s82_4_B','3u5m_27_C','3vot_3_A','4ecg_5_A','4hap_5_A','4jcl_4_A','4ktp_1_B','4mhv_1_A','4oy7_5_FH','4xb3_3_A','5jcw_2_B'],
    5:  ['1g5c_4_B','2ozb_2_E','1mjo_1_B','1t9h_6_A','2awy_3_AB','2b50_2_A','2h9d_1_CD','2x49_3_A','3b40_1_A','3cq8_7_A','3g0k_1_A','3qhq_2_B','3qta_2_A','4a15_2_A','4c2a_2_A','4ex8_1_A','4n17_2_A','4yl9_1_D','5klb_2_ABD'],
    6:  ['1jn9_3_A','3mhr_3_A','1o5k_1_A','1scf_2_C','2hn2_1_AB','2j45_2_AB','2jcg_1_A','2obl_1_A','2vyo_1_A','3aua_2_B','3faw_2_A','3hjb_11_D','3o83_5_A','3ohe_2_B','3tb3_1_B','5i29_1_A','5jix_4_A','5l7p_1_A','5vkw_2_B'],
    7:  ['1q39_4_A','2xco_1_A','2bue_1_A','2iim_1_A','2p5r_1_A','3hdb_1_A','3ipr_3_F','3ktb_1_D','3np5_1_D','3u4l_2_AP','3uf5_1_A','3wl4_4_B','4apb_1_A','4dxw_1_A','5gmf_2_D'],
    8:  ['1w3b_1_A','4lw9_48_Q','2bh1_1_A','2wq5_1_A','3bvc_6_A','3gri_6_A','3gv5_1_D','3lez_4_A','3r40_1_A','3tqk_2_A','3whn_1_A','4hhr_5_A','4j3h_3_A','4uix_1_B','6oz7_2_H'],
    9:  ['1y10_7_CD','5ue1_3_A','2b01_2_A','2buk_3_A','3a51_7_A','3jq1_8_B','3lez_6_A','3mdo_4_B','3mzo_4_A','3upg_2_A','4jbe_1_A','4r36_1_AB','5az1_1_A','5ts2_1_AD','5uxa_1_A','6gsz_9_A'],
    10: ['2c9m_3_A','5y8x_3_A','2xts_12_AB','2zal_4_AB','3jq1_9_B','3lv4_2_A','3o83_4_B','4kn9_15_T','4mbo_2_A','4q6b_3_A','5cyb_1_A'],
    11: ['2j45_1_A','4m0k_3_D','3iqt_1_A','4b5w_2_C','4eqb_4_A','4pih_3_B','5a8k_5_A','5bqt_1_BD','5nsf_1_A','6uff_12_D'],
}

print(f"{'Cluster':>8} {'Pairs':>7} {'d=1.0':>7} {'d=1.0%':>8} {'<0.5':>6} {'<0.5% (all)':>12} {'<0.5% (excl 1.0)':>18}")
print("-" * 80)

total_pairs = total_missing = total_below = 0

for cid, members in clusters.items():
    valid = [m for m in members if m in df.index]
    pairs_total = pairs_missing = pairs_below = 0
    for i in range(len(valid)):
        for j in range(i+1, len(valid)):
            d = df.loc[valid[i], valid[j]]
            pairs_total += 1
            if d == 1.0:
                pairs_missing += 1
            if d < 0.5:
                pairs_below += 1
    pct_all    = 100 * pairs_below / pairs_total if pairs_total else 0
    denom_excl = pairs_total - pairs_missing
    pct_excl   = 100 * pairs_below / denom_excl if denom_excl else 0
    print(f"{cid:>8} {pairs_total:>7} {pairs_missing:>7} {100*pairs_missing/pairs_total:>7.1f}% {pairs_below:>6} {pct_all:>11.1f}% {pct_excl:>17.1f}%")
    total_pairs   += pairs_total
    total_missing += pairs_missing
    total_below   += pairs_below

print("-" * 80)
pct_all_total  = 100 * total_below / total_pairs
pct_excl_total = 100 * total_below / (total_pairs - total_missing)
print(f"{'TOTAL':>8} {total_pairs:>7} {total_missing:>7} {100*total_missing/total_pairs:>7.1f}% {total_below:>6} {pct_all_total:>11.1f}% {pct_excl_total:>17.1f}%")
print(f"\n=> pct_below_0.5 incluindo d=1.0:  {pct_all_total:.2f}%  (threshold_analysis.py)")
print(f"=> pct_below_0.5 excluindo d=1.0:  {pct_excl_total:.2f}%  (possível origem da Table 4)")
