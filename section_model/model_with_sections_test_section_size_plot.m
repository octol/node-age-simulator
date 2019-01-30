clear all
load nodes01.dat
for ii = 1:length(number_of_sections)
    section_size = sum(nodes{ii}.active, 2);

    malicious_per_section{ii} = sum(nodes{ii}.malicious.*nodes{ii}.active, 2);
    malicious_fraction = malicious_per_section{ii} ./ section_size;
    malicious_fraction_mean(ii) = mean(malicious_fraction);
    malicious_fraction_std(ii) = std(malicious_fraction);
    malicious_fraction_max(ii) = max(malicious_fraction);

    malicious_elder_per_section{ii} = sum(nodes{ii}.malicious.*nodes{ii}.elder.*nodes{ii}.active, 2);
    malicious_elder_fraction = malicious_elder_per_section{ii} ./ section_size;
    malicious_elder_fraction_mean(ii) = mean(malicious_elder_fraction);
    malicious_elder_fraction_std(ii) = std(malicious_elder_fraction);
    malicious_elder_fraction_max(ii) = max(malicious_elder_fraction);
end

figure(1); clf;
hold on
plot(min_section_size, malicious_fraction_mean, 'LineWidth', 2);
plot(min_section_size, malicious_fraction_mean + malicious_fraction_std, 'LineWidth', 2);
plot(min_section_size, malicious_fraction_max, 'LineWidth', 2);
hold off
xlabel('Section size')
ylabel('Fraction')
legend({'\sigma', '\mu', 'max'})
title('Malicious node / section');

figure(2); clf;
hold on
plot(min_section_size, malicious_elder_fraction_mean, 'LineWidth', 2);
plot(min_section_size, malicious_elder_fraction_mean + malicious_elder_fraction_std, 'LineWidth', 2);
plot(min_section_size, malicious_elder_fraction_max, 'LineWidth', 2);
hold off
xlabel('Section size')
ylabel('Fraction')
legend({'\sigma', '\mu', 'max'})
title('Malicious elder / section');

%figure(2); clf;
%hist(malicious_per_section{end},100);
