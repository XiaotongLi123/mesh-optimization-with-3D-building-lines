function [intersect_section]=intersectSection(section_1,section_2)
    intersect_section=[];
    section_1=sort(section_1);
    section_2=sort(section_2);
    if section_1(2,1)<section_2(1,1)||section_2(2,1)<section_1(1,1)
        return;
    else
        point_sort=sort([section_1;section_2]);
        intersect_section=[point_sort(2,1);point_sort(3,1)];
    end
end
