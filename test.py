from hijri_converter import convert

# Gregorian date to convert
gregorian_date = '2023-06-23'

# Convert Gregorian to Hijri
hijri_date = convert.Gregorian(gregorian_date).to_hijri()

print(hijri_date)