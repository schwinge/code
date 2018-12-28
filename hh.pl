=https://www.learn-perl.org/en/Conditional_Decisions
An array @family holds a list of family member names. The first hash %shoe_color contains favorite shoe color per person name. The second hash %shoe_size contains shoe size per person name.

Evaluate and print the favorite shoe color and shoe size per each family member. For shoe sizes 10 and above, add the word 'large' to the output line.

Output lines should be in the format: "Homer wears large brown shoes size 12".

Note: not all family members may be included in the hash variables, so you better conditionally check if they exist or not (using the exists operator). If a name does not exist, add the key/value pair into the hash variables - for show color add: black; for shoe size add 99.
=cut
@family = ('Homer', 'Moe', 'Maggie');
%shoe_color = ('Lisa' => 'red', 'Homer' => 'brown', 'Maggie' => 'pink', 'Marge' => 'blue', 'Bart' => 'yellow');
%shoe_size = ('Moe' => 9, 'Lisa' => 7, 'Homer' => 12, 'Bart' => 8, 'Maggie' => 4);
# start your code here
if (exists $shoe_color{$family[0]} && $shoe_size{$family[0]} >= 10){
	print "$family[0] wears large $shoe_color{$family[0]} shoes size $shoe_size{$family[0]}\n";
} elsif(exists $shoe_color{$family[0]} && exists $shoe_size{$family[0]}){
	print "$family[0] wears $shoe_color{$family[0]} shoes size $shoe_size{$family[0]}\n";
} elsif(exists $shoe_color{$family[0]} && !(exists $shoe_size{$family[0]})){
    $shoe_size{$family[0]} = 99;
    print "$family[0] wears $shoe_color{$family[0]} shoes size $shoe_size{$family[0]}\n"
} else{
	$shoe_color{$family[0]} = black;
	$shoe_size{$family[0]} = 99;
    print "$family[0] wears $shoe_color{$family[0]} shoes size $shoe_size{$family[0]}\n"
}
if (exists $shoe_color{$family[1]} && $shoe_size{$family[1]} >= 10){
	print "$family[1] wears large $shoe_color{$family[1]} shoes size $shoe_size{$family[1]}\n";
} elsif(exists $shoe_color{$family[1]} && exists $shoe_size{$family[1]}){
	print "$family[1] wears $shoe_color{$family[1]} shoes size $shoe_size{$family[1]}\n";
} elsif(exists $shoe_color{$family[1]} && !(exists $shoe_size{$family[1]})){
    $shoe_size{$family[1]} = 99;
    print "$family[1] wears $shoe_color{$family[1]} shoes size $shoe_size{$family[1]}\n"
} else{
	$shoe_color{$family[1]} = black;
    print "$family[1] wears $shoe_color{$family[1]} shoes size $shoe_size{$family[1]}\n"
}
if (exists $shoe_color{$family[2]} && $shoe_size{$family[2]} >= 10){
	print "$family[2] wears large $shoe_color{$family[2]} shoes size $shoe_size{$family[2]}\n";
} elsif(exists $shoe_color{$family[2]} && exists $shoe_size{$family[2]}){
	print "$family[2] wears $shoe_color{$family[2]} shoes size $shoe_size{$family[2]}\n";
} elsif(exists $shoe_color{$family[2]} && !(exists $shoe_size{$family[2]})){
    $shoe_size{$family[2]} = 99;
    print "$family[2] wears $shoe_color{$family[2]} shoes size $shoe_size{$family[2]}\n"
} else{
	$shoe_color{$family[2]} = black;
    print "$family[2] wears $shoe_color{$family[2]} shoes size $shoe_size{$family[2]}\n"
}
=Homer wears large brown shoes size 12
Moe wears black shoes size 9
Maggie wears large pink shoes size 4
ge for string
=cut

@family = ('Homer', 'Moe', 'Maggie');
%shoe_color = ('Lisa' => 'red', 'Homer' => 'brown', 'Maggie' => 'pink', 'Marge' => 'blue', 'Bart' => 'yellow');
%shoe_size = ('Moe' => 9, 'Lisa' => 7, 'Homer' => 12, 'Bart' => 8, 'Maggie' => 4);

$default_shoe_color = "black";
$default_shoe_size = 4;

$member = $family[0];
if (!exists $shoe_color{$member}) {
	$shoe_color{$member} = $default_shoe_color;
}
if (!exists $shoe_size{$member}) {
	$shoe_size{$member} = $default_shoe_size;
}
$is_large = ($shoe_size{$member} >= 10) ? " large " : " ";
print "$member wears$is_large$shoe_color{$member} shoes size $shoe_size{$member}\n";

$member = $family[1];
if (!exists $shoe_color{$member}) {
	$shoe_color{$member} = $default_shoe_color;
}
if (!exists $shoe_size{$member}) {
	$shoe_size{$member} = $default_shoe_size;
}
$is_large = ($shoe_size{$member} >= 10) ? " large " : " ";
print "$member wears$is_large$shoe_color{$member} shoes size $shoe_size{$member}\n";

$member = $family[2];
if (!exists $shoe_color{$member}) {
	$shoe_color{$member} = $default_shoe_color;
}
if (!exists $shoe_size{$member}) {
	$shoe_size{$member} = $default_shoe_size;
}
$is_large = ($shoe_size{$member} >= 10) ? " large " : " ";
print "$member wears$is_large$shoe_color{$member} shoes size $shoe_size{$member}\n";