number_of_seconds_per_day = 86400.0


def seconds_to_days(duration_in_seconds):
    return duration_in_seconds / number_of_seconds_per_day


def days_to_seconds(duration_in_days):
    return duration_in_days * number_of_seconds_per_day


def per_day_to_per_second(rate_per_day):
    return rate_per_day / number_of_seconds_per_day


def per_second_to_per_day(rate_per_second):
    return rate_per_second * number_of_seconds_per_day


def Pa_to_kPa(stress_in_Pa):
    return stress_in_Pa / 1000.0


def unit_to_k_unit(unit_value):
    return unit_value / 1000.0


def fraction_to_percentage(fraction_value):
    return fraction_value * 100.0
