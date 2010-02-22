-- Autogenerated on Fri May  8 09:42:00 2009 by mkaudit.pl

create table audit.feature_cvtermprop (
    feature_cvtermprop_id integer not null
  , value text
  , type_id integer not null
  , feature_cvterm_id integer not null
  , rank integer not null
) inherits (audit.audit);

create or replace function audit.audit_feature_cvtermprop_insert_proc()
returns trigger
as $$
BEGIN
  raise exception 'Cannot insert directly into audit.feature_cvtermprop. Use one of the child tables.';
END;
$$ language plpgsql;
create trigger feature_cvtermprop_insert_tr before insert on audit.feature_cvtermprop
    for each statement execute procedure audit.audit_feature_cvtermprop_insert_proc();
grant select on audit.feature_cvtermprop to chado_ro_role;
grant select, insert on audit.feature_cvtermprop to chado_rw_role;
grant execute on function audit.audit_feature_cvtermprop_insert_proc() to chado_rw_role;


create table audit.feature_cvtermprop_insert (
    constraint feature_cvtermprop_insert_ck check (type = 'INSERT')
) inherits (audit.feature_cvtermprop);
alter table audit.feature_cvtermprop_insert alter type set default 'INSERT';
grant select on audit.feature_cvtermprop_insert to chado_ro_role;
grant select, insert on audit.feature_cvtermprop_insert to chado_rw_role;

create or replace function audit.public_feature_cvtermprop_insert_proc()
returns trigger
as $$
BEGIN
  insert into audit.feature_cvtermprop_insert (
      feature_cvtermprop_id, feature_cvterm_id, type_id, value, rank
  ) values (
      new.feature_cvtermprop_id, new.feature_cvterm_id, new.type_id, new.value, new.rank
  );
  return new;
END;
$$ language plpgsql;
create trigger feature_cvtermprop_audit_insert_tr after insert on public.feature_cvtermprop
    for each row execute procedure audit.public_feature_cvtermprop_insert_proc();
grant execute on function audit.public_feature_cvtermprop_insert_proc() to chado_rw_role;


create table audit.feature_cvtermprop_update (
    constraint feature_cvtermprop_update_ck check (type = 'UPDATE')
  , old_type_id integer not null
  , old_rank integer not null
  , old_value text
  , old_feature_cvterm_id integer not null
) inherits (audit.feature_cvtermprop);
alter table audit.feature_cvtermprop_update alter type set default 'UPDATE';
grant select on audit.feature_cvtermprop_update to chado_ro_role;
grant select, insert on audit.feature_cvtermprop_update to chado_rw_role;

create or replace function audit.public_feature_cvtermprop_update_proc()
returns trigger
as $$
BEGIN
  if old.feature_cvtermprop_id <> new.feature_cvtermprop_id or old.feature_cvtermprop_id is null <> new.feature_cvtermprop_id is null then
    raise exception 'If you want to change feature_cvtermprop.feature_cvtermprop_id (do you really?) then disable the audit trigger feature_cvtermprop_audit_update_tr';
  end if;
  insert into audit.feature_cvtermprop_update (
      feature_cvtermprop_id, feature_cvterm_id, type_id, value, rank,
      old_feature_cvterm_id, old_type_id, old_value, old_rank
   ) values (
       new.feature_cvtermprop_id, new.feature_cvterm_id, new.type_id, new.value, new.rank,
       old.feature_cvterm_id, old.type_id, old.value, old.rank
   );
  return new;
END;
$$ language plpgsql;
create trigger feature_cvtermprop_audit_update_tr after update on public.feature_cvtermprop
    for each row execute procedure audit.public_feature_cvtermprop_update_proc();
grant execute on function audit.public_feature_cvtermprop_update_proc() to chado_rw_role;


create table audit.feature_cvtermprop_delete (
    constraint feature_cvtermprop_delete_ck check (type = 'DELETE')
) inherits (audit.feature_cvtermprop);
alter table audit.feature_cvtermprop_delete alter type set default 'DELETE';
grant select on audit.feature_cvtermprop_delete to chado_ro_role;
grant select, insert on audit.feature_cvtermprop_delete to chado_rw_role;

create or replace function audit.public_feature_cvtermprop_delete_proc()
returns trigger
as $$
BEGIN
  insert into audit.feature_cvtermprop_delete (
      feature_cvtermprop_id, feature_cvterm_id, type_id, value, rank
  ) values (
      old.feature_cvtermprop_id, old.feature_cvterm_id, old.type_id, old.value, old.rank
  );
  return old;
END;
$$ language plpgsql;
create trigger feature_cvtermprop_audit_delete_tr after delete on public.feature_cvtermprop
    for each row execute procedure audit.public_feature_cvtermprop_delete_proc();
grant execute on function audit.public_feature_cvtermprop_delete_proc() to chado_rw_role;