args = commandArgs(trailingOnly = T)
dat = read.table(file ("stdin"), sep = "\t", row.names = 1, 
                 header = F, stringsAsFactors = F, 
                 comment.char = "%")
# Here all  valid states keep state IDs - 0, 1, and 2. It discards state `3` attributed to footprints on edges.
 
all_valid_states = data.frame(table(unlist(lapply (row.names(dat), 
                                           function (x) {as.integer(unlist(strsplit(x, split= "#"))[-1]) })))) 
names(all_valid_states) = c("state", "count")
state_info = c("2" = "N", "1" = "T", "0" = "D")
all_valid_states$state_label = unlist(lapply (all_valid_states$state, 
                                      function (x){state_info[as.character(x)]}))
all_valid_states$line_loc = rev(cumsum(rev(all_valid_states$count)))
all_valid_states$prev_point = c(all_valid_states$line_loc[-1], 0) 
all_valid_states$label_loc = (all_valid_states$line_loc + all_valid_states$prev_point)/2
all_valid_states$percentage = all_valid_states$count/sum (all_valid_states$count)
all_valid_states$total = sum (all_valid_states$count)
all_valid_states$annotation = args[1]
print (all_valid_states)
write.table(all_valid_states, paste0(args[1], ".tsv"), row.names = F, col.names = T, sep = "\t", quote = F)
total_molecules = nrow(dat)

jj = as.matrix (dat)
pdf(args[1], height = 6, width = 4.5)
par(mgp=c(1.5,0.25,0), cex = 0.75)
image(1:ncol(jj), 1:nrow(jj), t(jj),  axes = FALSE, useRaster = T, oldstyle = F, 
      col = c("-1" = "#bdbdbd", "0" = "#FFFFFF", "1"  = "Red" , "2" = "#2b8cbe"),
      xlim = c(-30, 301), xlab = args[2], ylab = "")
counter = 1
for (i in rev(seq (nrow(all_valid_states)) )){
    line_loc = all_valid_states[i, "line_loc"]
    label_loc = all_valid_states[i, "label_loc"]
    state_label = all_valid_states[i, "state_label"]
    cnt = all_valid_states[i, "count"]
    segments(-30, line_loc, 301, line_loc, lw = 0.5)
    if (counter == 1) {
        text(-15, label_loc, paste0 ("State ",  state_label, " ( n = ", cnt, ")"), srt = 90)
    }else{
        text(-15, label_loc, paste0 (state_label, " ( n = ", cnt, ")" ), srt = 90)
    }
    counter = counter + 1
}

# Draw vertical red line +/- 15 bases from `0`
segments(135, 0, 135, total_molecules, lw = 0.5, col = "Red") 
segments(165, 0, 165, total_molecules, lw = 0.5, col = "Red") 

axis(1, c(1, 50, 100, 150, 200, 250, 301) , c ("-150", "-100", "-50", "0", "50", "100", "150"), 
     srt = 2, cex.lab = 0.5, tck = -0.01, tick = F)
box(lwd=0.5)
dev.off()
