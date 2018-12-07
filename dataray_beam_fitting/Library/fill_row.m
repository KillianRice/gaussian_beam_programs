function out_table = fill_row(in_table, dist_cm, raw_data)


row_to_add = table;
row_to_add.Distance = dist_cm;
row_to_add.Gauss_u_mean_diam_um = mean(raw_data.G_2W);
row_to_add.Gauss_u_sig_diam_um = std(raw_data.G_2W);
row_to_add.Gauss_v_mean_diam_um = mean(raw_data.G_2W1);
row_to_add.Gauss_v_sig_diam_um = std(raw_data.G_2W1);
row_to_add.raw_data = {raw_data};

if isempty(in_table)
    out_table = row_to_add;
else
    in_table(height(in_table) + 1, :) = row_to_add;
    out_table = in_table;
end
end