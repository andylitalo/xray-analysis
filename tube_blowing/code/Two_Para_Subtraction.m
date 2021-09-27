function [alpha, beta] = Two_Para_Subtraction(img, bkg, q_low_min, q_low_max, q_high_min, q_high_max, q_matrix)
    
    %determine the low q and high q averages   
    counter1 = 1;
    counter2 = 1;
    im_low_q = 0;
    bkg_low_q = 0;
    im_high_q = 0;
    bkg_high_q = 0;
    
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            if q_matrix(i,j) >= q_low_min && q_matrix(i,j) <= q_low_max && img(i,j) > 0 && bkg(i,j) > 0
                im_low_q = im_low_q + img(i,j);
                bkg_low_q = bkg_low_q + bkg(i,j);
                counter1 = counter1 + 1;
            end
            if q_matrix(i,j) >= q_high_min && q_matrix(i,j) <= q_high_max && img(i,j) > 0 && bkg(i,j) > 0
                im_high_q = im_high_q + img(i,j);
                bkg_high_q = bkg_high_q + bkg(i,j);
                counter2 = counter2 + 1;
            end
        end
    end
    
    %compute low q and high q averages
    im_lowq_avg = im_low_q/counter1;
    im_highq_avg = im_high_q/counter2;
    bkg_lowq_avg = bkg_low_q/counter1;
    bkg_highq_avg = bkg_high_q/counter2;

    %compute scaling factors
    alpha = abs(im_highq_avg - im_lowq_avg)/abs(bkg_highq_avg - bkg_lowq_avg);
    beta = im_lowq_avg - alpha * bkg_lowq_avg;

end

